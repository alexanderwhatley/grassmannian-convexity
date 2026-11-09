// Check both event-separation witnesses mentioned in the report.
// This reuses the original checker's polynomial and derivative code; it is
// a reproducibility check, not a separate implementation of the model.
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wreturn-type"
#define main grassmann_original_main
#include "verify_upper_bound.cpp"
#undef main
#pragma GCC diagnostic pop

int main() {
    labels();
    auto events = symbolic_events(); // checks the witness in the main source
    vector<Poly> v;
    for (auto p : plabs) v.push_back(determinant(p[0], p[1], p[2]));
    for (auto f : flabs) v.push_back(feature({f[0],f[1]}, {f[2],f[3]}, {6,7}));
    v.push_back(feature({1,2}, {3,4}, {5,6}));
    v.push_back(feature({1,2}, {3,4}, {5,7}));
    assert(v.size()==52);
    array<array<int64_t,3>,7> M = {array<int64_t,3>{12,-5,1},
        {-7,4,-12}, {6,9,7}, {2,-11,-5}, {6,-11,-11},
        {-3,6,9}, {-6,12,-4}};
    array<int64_t,52> values{};
    for (int i=0;i<52;i++) { values[i]=evaluate(v[i],M); assert(values[i]); }
    bool active[7][52]={};
    for (auto e:events) active[e.move][e.target]=true;
    int checked=0;
    for (int j=1;j<=6;j++) {
        auto M1=M;
        for (int k=0;k<3;k++) M1[j-1][k]=add64(M1[j-1][k],M[j][k]);
        array<int64_t,52> d{};
        for (int i=0;i<52;i++) {
            d[i]=sub64(evaluate(v[i],M1),values[i]);
            assert((d[i]!=0)==active[j][i]);
        }
        for (int a=0;a<52;a++) if (active[j][a])
            for (int b=a+1;b<52;b++) if (active[j][b]) {
                assert(sub128(mul128(values[a],d[b]),mul128(values[b],d[a]))!=0);
                checked++;
            }
    }
    assert(checked==1070);
    cout << "PASS report witness: 52 features nonzero and " << checked
         << " same-generator resultants nonzero" << endl;
}
