//
// Created by Liangsheng Zhang on 9/18/15.
//

#ifndef MBL_V1_SORT_COMPARATOR_H
#define MBL_V1_SORT_COMPARATOR_H

/**
 ** This file contains various comparators that can be used for sorting
 **/

#include <utility>

template<typename T1, typename T2>
bool pair_first_less_comparator( const std::pair<T1,T2>& l1, const std::pair<T1,T2>& l2){
    return l1.first < l2.first;
}

#endif //MBL_V1_SORT_COMPARATOR_H
