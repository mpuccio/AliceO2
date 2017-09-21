set(MODULE_NAME "ITSReconstruction")

O2_SETUP(NAME ${MODULE_NAME})

set(SRCS
    src/TrivialClustererTask.cxx
    src/ClustererTask.cxx
    src/CookedTrack.cxx
    src/CookedTrackerTask.cxx
    src/CA/Track.cxx    
    )
set(HEADERS
    include/${MODULE_NAME}/TrivialClustererTask.h
    include/${MODULE_NAME}/ClustererTask.h
    include/${MODULE_NAME}/CookedTrack.h
    include/${MODULE_NAME}/CookedTrackerTask.h
    include/${MODULE_NAME}/CA/Track.h
    )
set(NO_DICT_SRCS # sources not for the dictionary
    src/TrivialClusterer.cxx
    src/CookedTracker.cxx
    src/CA/Cell.cxx
    src/CA/Cluster.cxx
    src/CA/Event.cxx
    src/CA/IndexTable.cxx
    src/CA/Layer.cxx
    src/CA/MathUtils.cxx
    src/CA/PrimaryVertexContext.cxx
    src/CA/Road.cxx
    src/CA/Tracker.cxx
    src/CA/Tracklet.cxx
    src/CA/Vertexer.cxx
    src/CookedTracker.cxx
    )
set(NO_DICT_HEADERS # sources not for the dictionary
    include/${MODULE_NAME}/TrivialClusterer.h
    include/${MODULE_NAME}/CookedTracker.h
    include/${MODULE_NAME}/CA/Cell.h
    include/${MODULE_NAME}/CA/Cluster.h
    include/${MODULE_NAME}/CA/Event.h
    include/${MODULE_NAME}/CA/IndexTable.h
    include/${MODULE_NAME}/CA/Layer.h
    include/${MODULE_NAME}/CA/PrimaryVertexContext.h
    include/${MODULE_NAME}/CA/Road.h
    include/${MODULE_NAME}/CA/Tracker.h
    include/${MODULE_NAME}/CA/Tracklet.h
    include/${MODULE_NAME}/CA/Vertexer.h
    )
Set(LINKDEF src/ITSReconstructionLinkDef.h)
Set(LIBRARY_NAME ${MODULE_NAME})
Set(BUCKET_NAME its_reconstruction_bucket)
O2_GENERATE_LIBRARY()

