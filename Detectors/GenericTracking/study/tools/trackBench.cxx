// Copyright 2019-2020 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file trackBench.cxx
/// \brief o2-gtrk-track-bench: the tracker on a tracking input the sink
///        dumped, without DPL.
///
/// Reads what the sink wrote with `--dump-tracking-input`, builds the tracker
/// as the sink builds it, and runs it --repeat times on the host or, with
/// `--backend cuda|hip`, on a GPU through the backend's device tracker
/// (`--host-reference` then compares the last run's tracks with the host's).
/// Each run prints its digests, its score against the MC truth and the profile
/// of every phase; `--json` appends one line per run for the scripts, and the
/// studies the options ask for watch the last run.

#include "TrackBench.h"

#include "GenericTracking/Alice/TRK.h"

#include "CCDB/BasicCCDBManager.h"
#include "CommonUtils/NameConf.h"
#include "DataFormatsParameters/GRPMagField.h"

#include <algorithm>
#include <cctype>
#include <format>
#include <iostream>
#include <stdexcept>
#include <string>

namespace o2::gtrk::study
{

/// Load the field of the run's GRPMagField (`GLO/Config/GRPMagField` from the
/// CCDB) into o2::base::Propagator, as a workflow's GRPGeomHelper would.
/// Returns the nominal Bz.
double loadField(const std::string& url)
{
  const auto* grp = o2::ccdb::BasicCCDBManager::instance().get<o2::parameters::GRPMagField>("GLO/Config/GRPMagField");
  if (grp == nullptr) {
    throw std::runtime_error("no GLO/Config/GRPMagField at " + url);
  }
  o2::base::Propagator::initFieldFromGRP(grp);
  return workflow::nominalBzKiloGauss();
}

int run(int argc, char** argv)
{
  Options options;
  bpo::options_description description("Track a dumped tracking input");
  description.add_options()                                                                                                             //
    ("help,h", "this message")                                                                                                          //
    ("detector,d", bpo::value(&options.detector)->required(), "its, mft or trk: the ITS, the MFT or ALICE 3's barrel on O2's clusters") //
    ("ccdb-url", bpo::value(&options.ccdbUrl)->default_value(o2::base::NameConf::getCCDBServer()),
     "the conditions' CCDB, e.g. file://<sample>/ccdb for its snapshots: the aligned geometry (GLO/Config/GeometryAligned), the field (GLO/Config/GRPMagField, not "
     "for trk) and the material LUT (GLO/Param/MatLUT, MFT/Param/MatLUT, TRK/Param/MatLUT), the sink's conditions") //
    ("input,i", bpo::value(&options.input)->default_value("tracking-input.bin"),
     "the sink's --dump-tracking-input file")                                                                                                          //
    ("truth", bpo::value(&options.truth)->default_value(""), "its MC truth; <input>.truth if empty and that exists")                                   //
    ("configKeyValues", bpo::value(&options.keyValues)->default_value(""), "GTRKReco<DET> keys to change on top of the parameters, as the workflow's") //
    ("config", bpo::value(&options.config)->default_value(""), "the parameters it ran with; <input>.ini if empty")                                     //
    ("threads,t", bpo::value(&options.threads)->default_value(-1), "the task arena's threads; the configuration's GTRKReco<DET>.threads if negative")  //
    ("repeat,n", bpo::value(&options.repeat)->default_value(3U), "runs; the first is cold, the others warm")                                           //
    ("tree", bpo::value(&options.tree)->default_value(true), "print each run's phase tree")                                                            //
    ("iteration-trees", bpo::bool_switch(&options.iterationTrees), "print each iteration's phase tree too")                                            //
    ("no-profile", bpo::bool_switch(&options.noProfile), "run without the profile, to measure what it costs")                                          //
    ("dram", bpo::bool_switch(&options.dram), "count each phase's DRAM traffic (AMD memory-controller counters, system-wide: a quiet machine)")        //
    ("digest-stages", bpo::bool_switch(&options.digestStages), "print every slice's stage-output digests (an observer; it adds to the time)")          //
    ("findable-min-pt", bpo::value(&options.findableMinPtGeV)->default_value(0.1), "the score's findable floor [GeV]")                                 //
    ("findable-all-layers", bpo::value(&options.findableAllLayers)->default_value(0U),
     "the layers the score's \"on every layer\" classes ask, a bit per layer; 0 for all") //
    ("seed-coverage", bpo::value(&options.seedCoverage)->default_value(""),
     "layer masks, comma separated: the findable with a measurement on every layer of each, and on all but one") //
    ("seeding-strategy", bpo::value(&options.seedingStrategy)->default_value(""),
     "instead of the runs: each iteration's graph alone, on the rows the earlier ones' seeded particles left, its coverage and cost; one JSON line per pass "
     "here") //
    ("seeding-steps", bpo::value(&options.seedingSteps)->default_value("1.0,0.7,0.5,0.4,0.3,0.25,0.2,0.15,0.13,0.12,0.11,0.1"),
     "with --seeding-strategy: the pT steps its coverage is counted at [GeV]")                                                                  //
    ("seeding-no-removal", bpo::bool_switch(&options.seedingNoRemoval), "with --seeding-strategy: hold no rows, every pass on the whole frame") //
    ("seeding-shadow", bpo::value(&options.seedingShadow)->default_value(""),
     "with --seeding-strategy: these passes (comma separated) hold nothing, and the pass after each is compared with it on the same rows") //
    ("seeding-fit", bpo::bool_switch(&options.seedingFit),
     "with --seeding-strategy: a particle is seeded by a pure record of the road fit, not by a pure road") //
    ("seeding-real", bpo::value(&options.seedingReal)->default_value(""),
     "the real run's iterations scored as --seeding-strategy's passes (good tracks seed, the used mask holds), one JSON line per iteration here") //
    ("save-state", bpo::value(&options.saveState)->default_value(""),
     "run only the first --state-iterations iterations and write what they leave (tracks, used mask, summaries) here, for --resume-state") //
    ("state-iterations", bpo::value(&options.stateIterations)->default_value(0U), "with --save-state: the iterations it runs")             //
    ("hold-clusters", bpo::value(&options.holdClusters)->default_value(""),
     "a file of cluster numbers, one per line, whose rows the run takes before its first iteration: the host tracker, the observers on "
     "(the vertex study's --vertex-replay <prefix>-owned.txt)") //
    ("resume-state", bpo::value(&options.resumeState)->default_value(""),
     "go on from this --save-state checkpoint: its iterations must be the configuration's first; the host tracker, no observers; with "
     "--seeding-strategy the real run before the study's passes") //
    ("seeding-particles", bpo::value(&options.seedingParticles)->default_value(""),
     "with --seeding-strategy or --seeding-real: per findable particle (--seeding-population), the passes that seeded it, as CSV here") //
    ("seeding-population", bpo::value(&options.seedingPopulation)->default_value("primary"),
     "with --seeding-strategy: the findable particles it scores: primary (physical primaries), decay (decay daughters), v0 (those of a K0S or a "
     "hyperon)") //
    ("step-anatomy", bpo::value(&options.stepAnatomy)->default_value(""),
     "with --complete-anatomy: a CSV of every short track's refused step and of the complete tracks' inner steps, replayed") //
    ("hits", bpo::value(&options.hits)->default_value(""),
     "with --complete-anatomy: the MC hits (sgn_HitsTRK.root): the anatomy's path order, true kinks and every study against the MC")          //
    ("ft3-hits", bpo::value(&options.ft3Hits)->default_value("sgn_HitsFT3.root"), "TRKFT3: FT3's MC hits beside --hits (the disc crossings)") //
    ("kine", bpo::value(&options.kine)->default_value(""),
     "--complete-anatomy with --hits: the MC kinematics, a <prefix>_Kine.root (sgn_Kine.root), for the clutter test of the follow's steps") //
    ("truth-walk", bpo::value(&options.truthWalk)->default_value(""),
     "with --complete-anatomy and --hits: a CSV of every step of the truth walks, the step's prediction at the particle's own next row against the cluster and "
     "the MC crossing") //
    ("truth-walk-iteration", bpo::value(&options.truthWalkIteration)->default_value(-1),
     "with --complete-anatomy and --hits: the iteration the truth walk runs on (-1: the first whose layers hold the class's mask)") //
    ("inward-history", bpo::value(&options.inwardHistory)->default_value(""),
     "with --complete-anatomy and --hits: a CSV of the inward refit's states, the failing passes' and a control's, against the MC truth") //
    ("end-decisions", bpo::value(&options.endDecisions)->default_value(""),
     "with --complete-anatomy: a CSV of every follow walk's end -- how it ended, the unconfirmed drop's verdict, its last and dropped rows against the truth, "
     "whether the last row stayed on the track") //
    ("end-edits", bpo::value(&options.endEdits)->default_value(""),
     "with --end-decisions: a CSV of each track as it is and its truth-free end edits -- the last row refused, the unconfirmed drop's row put back -- "
     "refitted and floored") //
    ("anatomy-findable", bpo::bool_switch(&options.anatomyFindable),
     "with --complete-anatomy: only the particles findable by the QA's rule (a chain of 4 layers, --findable-min-pt)") //
    ("anatomy-beam-line", bpo::value(&options.anatomyBeamLine)->default_value(""),
     "with --complete-anatomy: a CSV of the beam-line test on the class's own roads, every iteration with the beam constraint, against the MC") //
    ("anatomy-follow-loss", bpo::value(&options.anatomyFollowLoss)->default_value(""),
     "with --complete-anatomy and --hits: a CSV of the class's selected records dropped or kept as a fake, against the MC momentum along the leg") //
    ("anatomy-pdg", bpo::value(&options.anatomyPdg)->default_value(0),
     "with --complete-anatomy: the class's species, |PDG code|; 0, every species") //
    ("breakdown", bpo::value(&options.breakdown)->default_value(""),
     "a JSON line appended: the last run's tracks by length (fakes, clones, complete), the findable particles of --breakdown-origin complete, incomplete "
     "and where the missing layers' rows went") //
    ("breakdown-particles", bpo::value(&options.breakdownParticles)->default_value(""),
     "a CSV appended: --breakdown's findable particles one by one, their layers, status, best track (layers, iteration) and missing layers by fate and "
     "place; a particle not found, every layer's fate") //
    ("breakdown-origin", bpo::value(&options.breakdownOrigin)->default_value("primary"),
     "--breakdown's class: primary (physical primaries) or decay (decay daughters)") //
    ("fate-out", bpo::value(&options.fateOut)->default_value(""),
     "with --complete-anatomy: a CSV of every class particle's fate after the last iteration, for a comparison between runs") //
    ("complete-anatomy", bpo::bool_switch(&options.completeAnatomy),
     "the findable on every layer of --findable-all-layers: found whole, short and why, not found and why") //
    ("follow-stops", bpo::bool_switch(&options.followStops),
     "the last run's follow ends of the correct roads against the truth, by direction and why")                                                      //
    ("funnel", bpo::bool_switch(&options.funnel), "the last run's stage funnel of the findable against the truth (observers; they add to the time)") //
    ("beam-width", bpo::value(&options.beamWidthCM)->default_value(0.),
     "the beam's transverse width the windows and the vertexer take [cm], instead of the dump's beam line's; 0 for the dump's") //
    ("vertex-study", bpo::value(&options.vertexStudy)->default_value(-1),
     "the last run's iteration whose triplets the vertex study takes, its vertexer against the truth and its replay; -1 for none") //
    ("vertices-only", bpo::bool_switch(&options.verticesOnly),
     "the vertexing pass alone (GTRKReco<DET>.vertexer), no tracking: every run timed, the first's vertices against the MC collisions") //
    ("score-sources", bpo::value(&options.scoreSources)->multitoken(),
     "the score again over the particles of these MC sources alone (the collisions added to a mixed timeframe), against every track") //
    ("vertex-list", bpo::value(&options.vertexList)->default_value(""),
     "the vertex study's vertexer vertices here, one CSV line each with its collision and makeup, for tools/vertexMakeup.py") //
    ("vertex-prior", bpo::value(&options.vertexPrior)->default_value(""),
     "an earlier pass's --vertex-replay <prefix>-vertices.csv: the replay's prior-* variants start each frame from its stage vertices") //
    ("vertex-replay", bpo::value(&options.vertexReplay)->default_value(""),
     "the vertex study's replay of the vertexer and its ownership variants: <prefix>-vertices.csv, <prefix>-collisions.csv, "
     "<prefix>-dropped.csv and <prefix>-owned.txt, for tools/vertexReplay.py")                                                                             //
    ("fit-inputs", bpo::value(&options.fitInputs)->default_value(""), "write the last run's triplet and neighbour fit inputs here, for o2-gtrk-fit-bench") //
    ("trace", bpo::value(&options.trace)->default_value(""), "write the last run's zones and chunks here, as a Chrome trace")                              //
    ("json", bpo::value(&options.json)->default_value(""), "append one JSON line per run here")                                                            //
    ("tracks-out", bpo::value(&options.tracksOut)->default_value(""), "write the last run's tracks here, one line each: iteration, clusters, chi2")        //
    ("follow-steps", bpo::value(&options.followSteps)->default_value(""),
     "write the follow's steps on correct roads against the truth here, replayed after each iteration") //
    ("fakes-out", bpo::value(&options.fakesOut)->default_value(""),
     "write the last run's tracks against the truth here: where each one's wrong rows are, road, follow or the earlier iterations' follow") //
    ("clones-out", bpo::value(&options.clonesOut)->default_value(""),
     "write the last run's findable physical primaries, their good tracks and the cloned ones' rows here") //
    ("rows-out", bpo::value(&options.rowsOut)->default_value(""),
     "write every row of the last run's frame here with its first particle's origin, kinematics and production point") //
    ("clone-joins", bpo::value(&options.cloneJoins)->default_value(""),
     "write the last run's clones joined along the truth here: each inner good track's outward state stepped through the next good track's rows, the true "
     "row's chi2 and the rows competing with it")                                                                                               //
    ("backend", bpo::value(&options.backend)->default_value("cpu"), "cpu, or cuda or hip: run each planned timeframe on that GPU backend")      //
    ("gpu-device", bpo::value(&options.gpuDevice)->default_value(0), "the GPU backend's device index")                                          //
    ("gpu-arena-mb", bpo::value(&options.gpuArenaMB)->default_value(0U), "the device arena [MB], in O2's GPU pool; a GPU backend needs a size") //
    ("host-reference", bpo::bool_switch(&options.hostReference), "after the GPU runs, track on the host once and compare the tracks")           //
    ("rebin", bpo::bool_switch(&options.rebin), "rebin the rows even when they are binned by this grid already, which must change nothing");

  bpo::variables_map vm;
  try {
    bpo::store(bpo::parse_command_line(argc, argv, description), vm);
    bpo::notify(vm);
  } catch (const bpo::error& error) {
    std::cerr << error.what() << "\n\n"
              << description << '\n';
    return 1;
  }
  if (vm.contains("help")) {
    std::cout << description << '\n';
    return 0;
  }

  // The trial's keys and --threads go on in track, after the dump's own grid
  // is known.
  // The detector's parameter: GTRKRecoITS, GTRKRecoMFT, GTRKRecoTRK.
  std::string detector = options.detector;
  std::ranges::transform(detector, detector.begin(), [](unsigned char c) { return static_cast<char>(std::toupper(c)); });
  loadConfiguration(options.config.empty() ? options.input + ".ini" : options.config, "GTRKReco" + detector);
  // The conditions' CCDB: TrackBench.h reads the geometry and the material LUT from it too.
  o2::ccdb::BasicCCDBManager::instance().setURL(options.ccdbUrl);
  double bz = 0.;
  if (options.detector == "trk") {
    // ALICE 3's field is not in its GRP: O2's ALICE 3 field, as the workflow's conditions take it.
    bz = alice::alice3BzKiloGauss();
    o2::base::Propagator::Instance(true)->setNominalBz(static_cast<float>(bz));
    std::cout << std::format("the field: nominal Bz {} kG, O2's ALICE 3 field\n", bz);
  } else {
    bz = loadField(options.ccdbUrl);
    std::cout << std::format("the field: nominal Bz {} kG, from {}\n", bz, options.ccdbUrl);
  }
  if (options.detector == "its") {
    return track<detectors::ITS>(options);
  }
  if (options.detector == "mft") {
    return track<detectors::MFT>(options);
  }
  if (options.detector == "trk") {
    return track<detectors::TRKFT3>(options);
  }
  throw std::runtime_error("no detector called '" + options.detector + "'; its, mft or trk");
}

} // namespace o2::gtrk::study

int main(int argc, char** argv)
{
  try {
    return o2::gtrk::study::run(argc, argv);
  } catch (const std::exception& error) {
    std::cerr << "error: " << error.what() << '\n';
    return 1;
  }
}
