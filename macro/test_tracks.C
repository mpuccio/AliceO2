
void test_tracks() {
  float par[5] = {1.f,2.f,0.4f,0.5f,5.f};
  float cov[15] = {0.f};
  AliceO2::Base::Track::TrackPar a(0.f,0.f,par);
  AliceO2::Base::Track::TrackParCov c(0.f,0.f,par,cov);
  printf("(%f,%f,%f,%f,%f,%f,%f)\n",a[0],a[1],a[2],a[3],a[4],a[5],a[6]);
  AliceO2::Base::Track::TrackPar& d = c;
  printf("(%f,%f,%f,%f,%f,%f,%f)\n",d[0],d[1],d[2],d[3],d[4],d[5],d[6]);
  AliceO2::Base::Track::TrackPar e0 = c;
  printf("(%f,%f,%f,%f,%f,%f,%f)\n",e0[0],e0[1],e0[2],e0[3],e0[4],e0[5],e0[6]);
  PropagateParamTo(a,1.f, 0.5f);
  printf("(%f,%f,%f,%f,%f,%f,%f)\n",a[0],a[1],a[2],a[3],a[4],a[5],a[6]);
  PropagateParamTo(c,1.f, 0.5f);
  printf("(%f,%f,%f,%f,%f,%f,%f)\n",c[0],c[1],c[2],c[3],c[4],c[5],c[6]);
  
  for (int i = 0; i < AliceO2::Base::Track::kTrackPSize; ++i) {
    if (fabs(a[i] - c[i]) > AliceO2::Base::Track::kAlmost0) return;
  }
  cout << "Macro finished succesfully." << endl;
}
