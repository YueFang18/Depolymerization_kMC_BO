// =================================================================
// 修复后的 src/depolymerization_explicit_sequence_record_number.cpp
// =================================================================

#include "depolymerization_explicit_sequence_record_number.h"

// generate_two_different_indices函数的完整实现
std::pair<int, int> generate_two_different_indices(int range)
{
    int r1 = static_cast<int>(my_rand(range));
    int r2;
    do {
        r2 = static_cast<int>(my_rand(range));
    } while (r1 == r2);
    return {r1, r2};
}

void depolymerization_explicit_sequence_record_number(int reaction_index,
                                                      double *num_monomer,
                                                      int &num_chain,
                                                      vector<int> &M1M1,
                                                      vector<int> &M1M2,
                                                      vector<int> &M2M1,
                                                      vector<int> &M2M2,
                                                      poly_chain &poly,
                                                      double *tot_dyad,
                                                      double *tot_triad,
                                                      double monomerYield) {

    int randomPolymerIndex, randomPolymerIndex2, fissionChain1, fissionChain2;
    std::pair<int, int> SelectedChains;
    int r1, r2;
    vector<int>::iterator p_chain;

    static std::mt19937_64 rng{ std::random_device{}() };

    switch (reaction_index) {
        /************** Fission at unsaturated chain end **************/
        case 0:              //P_i = R_i-1 + R_1.
            species[0] -= 1;
            species[3] += 2;
            randomPolymerIndex = static_cast<int>(my_rand(unsaturated_Polymers.num_M1.size()));

            fissionChain1 = unsaturated_Polymers.num_M1[randomPolymerIndex] - 1;
            Pr.push_back({fissionChain1, 0});
            Pr.push_back({1, 1});

            p_chain = unsaturated_Polymers.num_M1.begin() + randomPolymerIndex;
            unsaturated_Polymers.num_M1.erase(p_chain);
            break;

            /**************** Fission in a saturated chain *************/
        case 1:             //P_i+j = R_i + R_j
            species[1] -= 1;
            species[3] += 2;

            // 按链长加权选择聚合物
            {
                double total_length = std::accumulate(
                        saturated_Polymers.num_M1.begin(),
                        saturated_Polymers.num_M1.end(),
                        0.0
                );

                std::uniform_real_distribution<double> uni_dist(0.0, total_length);
                double r = uni_dist(rng);

                double cum = 0.0;
                size_t weightedIndex = 0;
                for (size_t i = 0; i < saturated_Polymers.num_M1.size(); ++i) {
                    cum += static_cast<double>(saturated_Polymers.num_M1[i]);
                    if (r <= cum) {
                        weightedIndex = i;
                        break;
                    }
                }
                randomPolymerIndex = static_cast<int>(weightedIndex);
            }

            if ((saturated_Polymers.num_M1[randomPolymerIndex]) == 1) {
                num_monomer[0] += 1;
                species[4] += 1;
            } else {
                do {
                    fissionChain1 = my_rand(saturated_Polymers.num_M1[randomPolymerIndex]);
                    fissionChain2 = saturated_Polymers.num_M1[randomPolymerIndex] - fissionChain1;
                } while (static_cast<int>(fissionChain1) == 0 || static_cast<int>(fissionChain2) == 0);

                Pr.push_back({fissionChain1,0});
                Pr.push_back({fissionChain2,0});
            }

            p_chain = saturated_Polymers.num_M1.begin() + randomPolymerIndex;
            saturated_Polymers.num_M1.erase(p_chain);
            break;

            /**************** Fission at head-head defect *************/
        case 2:             //P_i+j = R_i + R_j
            species[2] -= 1;
            species[3] += 2;
            randomPolymerIndex = static_cast<int>(my_rand(HH_Polymers.num_M1.size()));

            fissionChain1 = HH_Polymers.chain1[randomPolymerIndex];
            Pr.push_back({fissionChain1,0});
            fissionChain2 = HH_Polymers.chain2[randomPolymerIndex];
            Pr.push_back({fissionChain2,0});

            p_chain = HH_Polymers.chain1.begin() + randomPolymerIndex;
            HH_Polymers.chain1.erase(p_chain);
            p_chain = HH_Polymers.chain2.begin() + randomPolymerIndex;
            HH_Polymers.chain2.erase(p_chain);
            p_chain = HH_Polymers.num_M1.begin() + randomPolymerIndex;
            HH_Polymers.num_M1.erase(p_chain);
            break;

            /****************Depropagation by beta-scission *************/
        case 3:              //R_i = R_i-1 + M
            do {
                r1 = static_cast<int>(my_rand(Pr.size()));
            } while (Pr[r1][0] == 1);

            Pr[r1][0] -= 1;
            num_monomer[0] += 1;
            species[4] += 1;
            break;

            /****************ReCombination*************/
        case 4:              //Ri*+Rj*=P (H-H polymer)
            species[2] += 1;
            species[3] -= 2;
            std::tie(r1, r2) = generate_two_different_indices(Pr.size());

            HH_Polymers.chain1.push_back(Pr[r1][0]);
            HH_Polymers.chain2.push_back(Pr[r2][0]);
            HH_Polymers.num_M1.push_back(HH_Polymers.chain1.back() + HH_Polymers.chain2.back());

            if( r1 < r2 ) std::swap(r1, r2);
            Pr.erase(Pr.begin() + r1);
            Pr.erase(Pr.begin() + r2);
            break;

            /****************Disproportionation*************/
        case 5:              // Ri + Rj = Pi(saturated) + Pj(unsaturated)
            species[3] -= 2;
            std::tie(r1, r2) = generate_two_different_indices(Pr.size());

            if (Pr[r1][0] == 1 && Pr[r2][0] == 1) {
                num_monomer[0] += 2;
                species[4] += 2;
            }
            else if (Pr[r1][0] != 1 && Pr[r2][0] != 1) {
                if (rand() % 2 == 0) {
                    saturated_Polymers.num_M1.push_back(Pr[r1][0]);
                    unsaturated_Polymers.num_M1.push_back(Pr[r2][0]);
                } else {
                    saturated_Polymers.num_M1.push_back(Pr[r2][0]);
                    unsaturated_Polymers.num_M1.push_back(Pr[r1][0]);
                }
                species[1] += 1;
                species[0] += 1;
            }
            else {
                int polymer_index = (Pr[r1][0] != 1) ? r1 : r2;
                if (rand() % 2 == 0) {
                    unsaturated_Polymers.num_M1.push_back(Pr[polymer_index][0]);
                    species[0] += 1;
                } else {
                    saturated_Polymers.num_M1.push_back(Pr[polymer_index][0]);
                    species[1] += 1;
                }
                num_monomer[0] += 1;
                species[4] += 1;
            }

            if (r1 < r2) std::swap(r1, r2);
            Pr.erase(Pr.begin() + r1);
            Pr.erase(Pr.begin() + r2);
            break;
    }

    species[3] = Pr.size(); // update species[3]
    species[5] = 0;
    for (size_t i = 0; i < Pr.size(); ++i) {
        if (Pr[i][0] == 1) {
            species[5]++;
        }
    }
}