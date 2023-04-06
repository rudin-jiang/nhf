#include "hgp/hgp_hrr.hpp"
#include <vector>

namespace nhfInt {
namespace hgp {

void hgp_hrr_5_1_(std::vector<double> &a5b1, const std::vector<double> &hrrInp, double x, double y, double z)
{
    const double *a5b0 = hrrInp.data() +  0;
    const double *a6b0 = hrrInp.data() + 21;


    a5b1[ 0] = a6b0[ 0] + x * a5b0[ 0];    //[(500)(100)] = [(600)(000)] + ABx * [(500)(000)]
    a5b1[ 1] = a6b0[ 1] + y * a5b0[ 0];    //[(500)(010)] = [(510)(000)] + ABy * [(500)(000)]
    a5b1[ 2] = a6b0[ 2] + z * a5b0[ 0];    //[(500)(001)] = [(501)(000)] + ABz * [(500)(000)]
    a5b1[ 3] = a6b0[ 1] + x * a5b0[ 1];    //[(410)(100)] = [(510)(000)] + ABx * [(410)(000)]
    a5b1[ 4] = a6b0[ 3] + y * a5b0[ 1];    //[(410)(010)] = [(420)(000)] + ABy * [(410)(000)]
    a5b1[ 5] = a6b0[ 4] + z * a5b0[ 1];    //[(410)(001)] = [(411)(000)] + ABz * [(410)(000)]
    a5b1[ 6] = a6b0[ 2] + x * a5b0[ 2];    //[(401)(100)] = [(501)(000)] + ABx * [(401)(000)]
    a5b1[ 7] = a6b0[ 4] + y * a5b0[ 2];    //[(401)(010)] = [(411)(000)] + ABy * [(401)(000)]
    a5b1[ 8] = a6b0[ 5] + z * a5b0[ 2];    //[(401)(001)] = [(402)(000)] + ABz * [(401)(000)]
    a5b1[ 9] = a6b0[ 3] + x * a5b0[ 3];    //[(320)(100)] = [(420)(000)] + ABx * [(320)(000)]
    a5b1[10] = a6b0[ 6] + y * a5b0[ 3];    //[(320)(010)] = [(330)(000)] + ABy * [(320)(000)]
    a5b1[11] = a6b0[ 7] + z * a5b0[ 3];    //[(320)(001)] = [(321)(000)] + ABz * [(320)(000)]
    a5b1[12] = a6b0[ 4] + x * a5b0[ 4];    //[(311)(100)] = [(411)(000)] + ABx * [(311)(000)]
    a5b1[13] = a6b0[ 7] + y * a5b0[ 4];    //[(311)(010)] = [(321)(000)] + ABy * [(311)(000)]
    a5b1[14] = a6b0[ 8] + z * a5b0[ 4];    //[(311)(001)] = [(312)(000)] + ABz * [(311)(000)]
    a5b1[15] = a6b0[ 5] + x * a5b0[ 5];    //[(302)(100)] = [(402)(000)] + ABx * [(302)(000)]
    a5b1[16] = a6b0[ 8] + y * a5b0[ 5];    //[(302)(010)] = [(312)(000)] + ABy * [(302)(000)]
    a5b1[17] = a6b0[ 9] + z * a5b0[ 5];    //[(302)(001)] = [(303)(000)] + ABz * [(302)(000)]
    a5b1[18] = a6b0[ 6] + x * a5b0[ 6];    //[(230)(100)] = [(330)(000)] + ABx * [(230)(000)]
    a5b1[19] = a6b0[10] + y * a5b0[ 6];    //[(230)(010)] = [(240)(000)] + ABy * [(230)(000)]
    a5b1[20] = a6b0[11] + z * a5b0[ 6];    //[(230)(001)] = [(231)(000)] + ABz * [(230)(000)]
    a5b1[21] = a6b0[ 7] + x * a5b0[ 7];    //[(221)(100)] = [(321)(000)] + ABx * [(221)(000)]
    a5b1[22] = a6b0[11] + y * a5b0[ 7];    //[(221)(010)] = [(231)(000)] + ABy * [(221)(000)]
    a5b1[23] = a6b0[12] + z * a5b0[ 7];    //[(221)(001)] = [(222)(000)] + ABz * [(221)(000)]
    a5b1[24] = a6b0[ 8] + x * a5b0[ 8];    //[(212)(100)] = [(312)(000)] + ABx * [(212)(000)]
    a5b1[25] = a6b0[12] + y * a5b0[ 8];    //[(212)(010)] = [(222)(000)] + ABy * [(212)(000)]
    a5b1[26] = a6b0[13] + z * a5b0[ 8];    //[(212)(001)] = [(213)(000)] + ABz * [(212)(000)]
    a5b1[27] = a6b0[ 9] + x * a5b0[ 9];    //[(203)(100)] = [(303)(000)] + ABx * [(203)(000)]
    a5b1[28] = a6b0[13] + y * a5b0[ 9];    //[(203)(010)] = [(213)(000)] + ABy * [(203)(000)]
    a5b1[29] = a6b0[14] + z * a5b0[ 9];    //[(203)(001)] = [(204)(000)] + ABz * [(203)(000)]
    a5b1[30] = a6b0[10] + x * a5b0[10];    //[(140)(100)] = [(240)(000)] + ABx * [(140)(000)]
    a5b1[31] = a6b0[15] + y * a5b0[10];    //[(140)(010)] = [(150)(000)] + ABy * [(140)(000)]
    a5b1[32] = a6b0[16] + z * a5b0[10];    //[(140)(001)] = [(141)(000)] + ABz * [(140)(000)]
    a5b1[33] = a6b0[11] + x * a5b0[11];    //[(131)(100)] = [(231)(000)] + ABx * [(131)(000)]
    a5b1[34] = a6b0[16] + y * a5b0[11];    //[(131)(010)] = [(141)(000)] + ABy * [(131)(000)]
    a5b1[35] = a6b0[17] + z * a5b0[11];    //[(131)(001)] = [(132)(000)] + ABz * [(131)(000)]
    a5b1[36] = a6b0[12] + x * a5b0[12];    //[(122)(100)] = [(222)(000)] + ABx * [(122)(000)]
    a5b1[37] = a6b0[17] + y * a5b0[12];    //[(122)(010)] = [(132)(000)] + ABy * [(122)(000)]
    a5b1[38] = a6b0[18] + z * a5b0[12];    //[(122)(001)] = [(123)(000)] + ABz * [(122)(000)]
    a5b1[39] = a6b0[13] + x * a5b0[13];    //[(113)(100)] = [(213)(000)] + ABx * [(113)(000)]
    a5b1[40] = a6b0[18] + y * a5b0[13];    //[(113)(010)] = [(123)(000)] + ABy * [(113)(000)]
    a5b1[41] = a6b0[19] + z * a5b0[13];    //[(113)(001)] = [(114)(000)] + ABz * [(113)(000)]
    a5b1[42] = a6b0[14] + x * a5b0[14];    //[(104)(100)] = [(204)(000)] + ABx * [(104)(000)]
    a5b1[43] = a6b0[19] + y * a5b0[14];    //[(104)(010)] = [(114)(000)] + ABy * [(104)(000)]
    a5b1[44] = a6b0[20] + z * a5b0[14];    //[(104)(001)] = [(105)(000)] + ABz * [(104)(000)]
    a5b1[45] = a6b0[15] + x * a5b0[15];    //[(050)(100)] = [(150)(000)] + ABx * [(050)(000)]
    a5b1[46] = a6b0[21] + y * a5b0[15];    //[(050)(010)] = [(060)(000)] + ABy * [(050)(000)]
    a5b1[47] = a6b0[22] + z * a5b0[15];    //[(050)(001)] = [(051)(000)] + ABz * [(050)(000)]
    a5b1[48] = a6b0[16] + x * a5b0[16];    //[(041)(100)] = [(141)(000)] + ABx * [(041)(000)]
    a5b1[49] = a6b0[22] + y * a5b0[16];    //[(041)(010)] = [(051)(000)] + ABy * [(041)(000)]
    a5b1[50] = a6b0[23] + z * a5b0[16];    //[(041)(001)] = [(042)(000)] + ABz * [(041)(000)]
    a5b1[51] = a6b0[17] + x * a5b0[17];    //[(032)(100)] = [(132)(000)] + ABx * [(032)(000)]
    a5b1[52] = a6b0[23] + y * a5b0[17];    //[(032)(010)] = [(042)(000)] + ABy * [(032)(000)]
    a5b1[53] = a6b0[24] + z * a5b0[17];    //[(032)(001)] = [(033)(000)] + ABz * [(032)(000)]
    a5b1[54] = a6b0[18] + x * a5b0[18];    //[(023)(100)] = [(123)(000)] + ABx * [(023)(000)]
    a5b1[55] = a6b0[24] + y * a5b0[18];    //[(023)(010)] = [(033)(000)] + ABy * [(023)(000)]
    a5b1[56] = a6b0[25] + z * a5b0[18];    //[(023)(001)] = [(024)(000)] + ABz * [(023)(000)]
    a5b1[57] = a6b0[19] + x * a5b0[19];    //[(014)(100)] = [(114)(000)] + ABx * [(014)(000)]
    a5b1[58] = a6b0[25] + y * a5b0[19];    //[(014)(010)] = [(024)(000)] + ABy * [(014)(000)]
    a5b1[59] = a6b0[26] + z * a5b0[19];    //[(014)(001)] = [(015)(000)] + ABz * [(014)(000)]
    a5b1[60] = a6b0[20] + x * a5b0[20];    //[(005)(100)] = [(105)(000)] + ABx * [(005)(000)]
    a5b1[61] = a6b0[26] + y * a5b0[20];    //[(005)(010)] = [(015)(000)] + ABy * [(005)(000)]
    a5b1[62] = a6b0[27] + z * a5b0[20];    //[(005)(001)] = [(006)(000)] + ABz * [(005)(000)]
}

} // namespace (hgp)
} // namespace (nhfInt)
