//
// Created by 13206 on 2024/1/31.
//

#ifndef EFFPROPARECORDER_POLYMERFEATHURE_H
#define EFFPROPARECORDER_POLYMERFEATHURE_H
#include "../include/DataStructureAndConstant.h"
void countPolymerFeatures(const poly_chain& poly,
                          double& num_unsaturated,
                          double& num_saturated,
                          double& num_HHbond) {
    num_unsaturated = 0;
    num_saturated = 0;
    num_HHbond = 0;

    for (int i = 0; i < poly.chain_feature.size(); i++) {
        if (poly.chain_feature[i] == 0 or poly.chain_feature[i] == 00000) {
            num_unsaturated++;
        } else if (poly.chain_feature[i] == 1 or poly.chain_feature[i] == 11111) {
            num_saturated++;
        } else if (poly.chain_feature[i] == 2) {
            num_HHbond++;
        }
        else if (poly.chain_feature[i] == 00000) {
            num_unsaturated++;
            num_saturated++;
        }
        else if (poly.chain_feature[i] == 11111) {
            num_saturated++;
            num_saturated++;
        }
    }
}

#endif //EFFPROPARECORDER_POLYMERFEATHURE_H
