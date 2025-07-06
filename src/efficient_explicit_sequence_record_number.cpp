// =================================================================
// 修复后的 src/efficient_explicit_sequence_record_number.cpp
// =================================================================

#include "efficient_explicit_sequence_record_number.h"

// ChooseChain函数的完整实现
std::pair<int,int> ChooseChain(const std::vector<std::vector<int>>& Pr, double wp)
{
    static std::mt19937_64 rng{ std::random_device{}() };

    struct Bucket { std::vector<int> idx; };
    std::unordered_map<int, Bucket> bucket;
    bucket.reserve(Pr.size());

    for (int k = 0; k < static_cast<int>(Pr.size()); ++k)
        bucket[ Pr[k][0] ].idx.push_back(k);

    struct Cum { double acc; int Li, Lj; };
    std::vector<Cum> table;
    table.reserve(bucket.size() * bucket.size() / 2 + 8);

    double total = 0.0;

    for (auto it_i = bucket.begin(); it_i != bucket.end(); ++it_i)
    {
        const int Li = it_i->first;
        const int Ni = static_cast<int>(it_i->second.idx.size());

        if (Ni > 1)
        {
            const double kij = CalculateKt_ij(wp, Li, Li);
            total += kij * 0.5 * Ni * (Ni - 1);
            table.push_back({ total, Li, Li });
        }

        for (auto it_j = std::next(it_i); it_j != bucket.end(); ++it_j)
        {
            const int Lj = it_j->first;
            const int Nj = static_cast<int>(it_j->second.idx.size());
            if (Nj == 0) continue;

            const double kij = CalculateKt_ij(wp, Li, Lj);
            total += kij * Ni * Nj;
            table.push_back({ total, Li, Lj });
        }
    }

    if (total == 0.0) return { -1, -1 };

    std::uniform_real_distribution<double> Ureal(0.0, std::nextafter(total, std::numeric_limits<double>::lowest()));
    const double r = Ureal(rng);

    auto pick = std::lower_bound(
            table.begin(), table.end(), r,
            [](const Cum& c, double v) { return c.acc < v; });

    if (pick == table.end()) --pick;

    const int Li = pick->Li;
    const int Lj = pick->Lj;

    const auto& vec_i = bucket[Li].idx;
    const auto& vec_j = bucket[Lj].idx;

    int i, j;

    if (Li == Lj)
    {
        const int Ni = static_cast<int>(vec_i.size());
        std::uniform_int_distribution<long long> Ucomb(0LL, static_cast<long long>(Ni) * (Ni - 1) / 2 - 1);
        long long code = Ucomb(rng);

        int p = 0;
        while (code >= Ni - 1 - p)
        {
            code -= Ni - 1 - p;
            ++p;
        }
        const int q = p + 1 + static_cast<int>(code);

        i = vec_i[p];
        j = vec_i[q];
    }
    else
    {
        std::uniform_int_distribution<int> Ui(0, static_cast<int>(vec_i.size()) - 1);
        std::uniform_int_distribution<int> Uj(0, static_cast<int>(vec_j.size()) - 1);
        i = vec_i[Ui(rng)];
        j = vec_j[Uj(rng)];
    }

    if (i > j) std::swap(i, j);
    return { i, j };
}

void efficient_explicit_sequence_record_number(int reaction_index,
                                               double *num_monomer,
                                               double *species,
                                               int &num_chain,
                                               vector<int> &M1M1,
                                               vector<int> &M1M2,
                                               vector<int> &M2M1,
                                               vector<int> &M2M2,
                                               poly_chain &poly,
                                               double w_p) {
    int M1_counter = 0;
    std::pair<int, int> SelectedChains;
    static std::mt19937 rng{ std::random_device{}() };

    int r1, r2;
    vector<int>::iterator p_chain;

    switch (reaction_index) {
        /************** Initiator Decomposition **************/
        case 0:              //I2=2I*.
            species[0] -= 1;
            species[1] += 2;
            break;

            /**************** Chain Initiation*************/
        case 1 :             //I*+M1=I-M1*
            species[1] -= 1;
            species[2] -= 1;
            species[3] += 1;
            num_monomer[0] += 1;
            Pr.push_back({1,0}); //length=1;normal Pr;
            break;

            /**************** Chain Propagation*************/
        case 2 :             //Ri*+M1=Ri+1*
            species[2] -= 1;
            num_monomer[0] += 1;
            r1 = my_rand(Pr.size());
            Pr[r1][0]++;
            break;

            /****************ReCombination*************/
        case 3:              //Ri*+Rj*=P //H-H
            species[3] -= 2;
            species[4] += 1;

            SelectedChains = ChooseChain(Pr, w_p);
            r1 = SelectedChains.first;
            r2 = SelectedChains.second;
            randomlySwap(r1, r2);

            poly.chain1.push_back(Pr[r1][0]);
            poly.chain2.push_back(Pr[r2][0]);
            poly.num_M1.push_back(poly.chain1.back() + poly.chain2.back());
            poly.chain_feature.push_back(2); // head-head bond;
            poly.chain1_feature.push_back(Pr[r1][1]);
            poly.chain2_feature.push_back(Pr[r2][1]);

            if (Pr[r1][1]) {
                species[6] -= 1;
            }
            if (Pr[r2][1]) {
                species[6] -= 1;
            }

            if( r1 < r2 ) std::swap(r1, r2);
            Pr.erase(Pr.begin() +  r1 );
            Pr.erase(Pr.begin() +  r2);
            break;

            /****************Disproportionation*************/
        case 4:              //Ri*+Rj*=P+P
            species[3] -= 2;
            species[4] += 2;

            SelectedChains = ChooseChain(Pr, w_p);
            r1 = SelectedChains.first;
            r2 = SelectedChains.second;
            randomlySwap(r1, r2);

            poly.chain1.push_back(Pr[r1][0]);
            poly.chain2.push_back(0);
            poly.num_M1.push_back(poly.chain1.back() + poly.chain2.back());
            poly.chain_feature.push_back(0); // unsaturated chain end
            poly.chain1_feature.push_back(Pr[r1][1]);
            poly.chain2_feature.push_back(66666);

            poly.chain1.push_back(Pr[r2][0]);
            poly.chain2.push_back(0);
            poly.num_M1.push_back(poly.chain1.back() + poly.chain2.back());
            poly.chain_feature.push_back(1); // saturated chain end
            poly.chain1_feature.push_back(Pr[r2][1]);
            poly.chain2_feature.push_back(66666);

            if (Pr[r1][1]) {
                species[6] -= 1;
            }
            if (Pr[r2][1]) {
                species[6] -= 1;
            }

            if( r1 < r2 ) std::swap(r1, r2);
            Pr.erase(Pr.begin() + r1 );
            Pr.erase(Pr.begin() + r2 );
            break;

            /**************** Chain transfer to monomer*************/
        case 5 :             //Ri*+M1=Pi+M*
            species[2] -= 1;
            species[3] -= 1;
            species[4] += 1;
            species[5] += 1;
            num_monomer[0] += 1;
            r1 = my_rand(Pr.size());

            poly.chain1.push_back(Pr[r1][0]);
            poly.chain2.push_back(0);
            poly.num_M1.push_back(poly.chain1.back() + poly.chain2.back());

            if (Pr[r1][1] == 1) {
                poly.chain_feature.push_back(00000); // unsaturated chain by transfer
                species[6] -= 1;
            }
            if (Pr[r1][1] == 0) {
                poly.chain_feature.push_back(11111);  // saturated chain end by transfer
            }
            poly.chain1_feature.push_back(Pr[r1][1]);
            poly.chain2_feature.push_back(99999);

            Pr.erase(Pr.begin() + r1);
            break;

            /**************** Chain initiation by monomer radical*************/
        case 6 :             //M*+M=M-M2*
            species[5] -= 1;
            species[2] -= 1;
            species[3] += 1;
            species[6] += 1;
            num_monomer[0] += 1;
            Pr.push_back({2,1});
            break;
    }
}