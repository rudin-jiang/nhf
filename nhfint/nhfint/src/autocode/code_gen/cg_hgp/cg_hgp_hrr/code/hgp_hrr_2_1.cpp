#include "hgp/hgp_hrr.hpp"
#include <vector>

namespace nhfInt {
namespace hgp {

void hgp_hrr_2_1_(std::vector<double> &a2b1, const std::vector<double> &hrrInp, double x, double y, double z)
{
    const double *a2b0 = hrrInp.data() + 0;
    const double *a3b0 = hrrInp.data() + 6;


    a2b1[ 0] = a3b0[ 0] + x * a2b0[ 0];    //[(200)(100)] = [(300)(000)] + ABx * [(200)(000)]
    a2b1[ 1] = a3b0[ 1] + y * a2b0[ 0];    //[(200)(010)] = [(210)(000)] + ABy * [(200)(000)]
    a2b1[ 2] = a3b0[ 2] + z * a2b0[ 0];    //[(200)(001)] = [(201)(000)] + ABz * [(200)(000)]
    a2b1[ 3] = a3b0[ 1] + x * a2b0[ 1];    //[(110)(100)] = [(210)(000)] + ABx * [(110)(000)]
    a2b1[ 4] = a3b0[ 3] + y * a2b0[ 1];    //[(110)(010)] = [(120)(000)] + ABy * [(110)(000)]
    a2b1[ 5] = a3b0[ 4] + z * a2b0[ 1];    //[(110)(001)] = [(111)(000)] + ABz * [(110)(000)]
    a2b1[ 6] = a3b0[ 2] + x * a2b0[ 2];    //[(101)(100)] = [(201)(000)] + ABx * [(101)(000)]
    a2b1[ 7] = a3b0[ 4] + y * a2b0[ 2];    //[(101)(010)] = [(111)(000)] + ABy * [(101)(000)]
    a2b1[ 8] = a3b0[ 5] + z * a2b0[ 2];    //[(101)(001)] = [(102)(000)] + ABz * [(101)(000)]
    a2b1[ 9] = a3b0[ 3] + x * a2b0[ 3];    //[(020)(100)] = [(120)(000)] + ABx * [(020)(000)]
    a2b1[10] = a3b0[ 6] + y * a2b0[ 3];    //[(020)(010)] = [(030)(000)] + ABy * [(020)(000)]
    a2b1[11] = a3b0[ 7] + z * a2b0[ 3];    //[(020)(001)] = [(021)(000)] + ABz * [(020)(000)]
    a2b1[12] = a3b0[ 4] + x * a2b0[ 4];    //[(011)(100)] = [(111)(000)] + ABx * [(011)(000)]
    a2b1[13] = a3b0[ 7] + y * a2b0[ 4];    //[(011)(010)] = [(021)(000)] + ABy * [(011)(000)]
    a2b1[14] = a3b0[ 8] + z * a2b0[ 4];    //[(011)(001)] = [(012)(000)] + ABz * [(011)(000)]
    a2b1[15] = a3b0[ 5] + x * a2b0[ 5];    //[(002)(100)] = [(102)(000)] + ABx * [(002)(000)]
    a2b1[16] = a3b0[ 8] + y * a2b0[ 5];    //[(002)(010)] = [(012)(000)] + ABy * [(002)(000)]
    a2b1[17] = a3b0[ 9] + z * a2b0[ 5];    //[(002)(001)] = [(003)(000)] + ABz * [(002)(000)]
}

} // namespace (hgp)
} // namespace (nhfInt)
