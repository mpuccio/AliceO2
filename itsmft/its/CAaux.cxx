#include "CAaux.h"

using namespace AliceO2::ITS::CA;

Cell::Cell(int xx,int yy, int zz, int dd0, int dd1, float curv, float n[3])
: m1OverR(curv),
  md0(dd0),
  md1(dd1),
  mN(),
  mVector() {
  mVector.reserve(4);
  mVector.push_back(xx);
  mVector.push_back(yy);
  mVector.push_back(zz);
  mVector.push_back(1u);
  if(n) {
    mN[0] = n[0];
    mN[1] = n[1];
    mN[2] = n[2];
  }
}

bool Cell::Combine(Cell &neigh, int idd) {
  // From outside inward
  if (this->y() == neigh.z() && this->x() == neigh.y()) { // Cells sharing two points
    mVector.push_back(idd);
    if (neigh.GetLevel() + 1 > GetLevel()) {
      SetLevel(neigh.GetLevel() + 1u);
    }
    return true;
  }
  return false;
}

