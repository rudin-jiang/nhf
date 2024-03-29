#include "hgp/hgp_hrr.hpp"
#include <vector>

namespace nhfInt {
namespace hgp {

void hgp_hrr_3_2_(std::vector<double> &a3b2, const std::vector<double> &hrrInp, double x, double y, double z)
{
    const double *a3b0 = hrrInp.data() +  0;
    const double *a4b0 = hrrInp.data() + 10;
    const double *a5b0 = hrrInp.data() + 25;

    double a3b1[30];
    double a4b1[45];

    a3b1[ 0] = a4b0[ 0] + x * a3b0[ 0];    //[(300)(100)] = [(400)(000)] + ABx * [(300)(000)]
    a3b1[ 1] = a4b0[ 1] + y * a3b0[ 0];    //[(300)(010)] = [(310)(000)] + ABy * [(300)(000)]
    a3b1[ 2] = a4b0[ 2] + z * a3b0[ 0];    //[(300)(001)] = [(301)(000)] + ABz * [(300)(000)]
    a3b1[ 3] = a4b0[ 1] + x * a3b0[ 1];    //[(210)(100)] = [(310)(000)] + ABx * [(210)(000)]
    a3b1[ 4] = a4b0[ 3] + y * a3b0[ 1];    //[(210)(010)] = [(220)(000)] + ABy * [(210)(000)]
    a3b1[ 5] = a4b0[ 4] + z * a3b0[ 1];    //[(210)(001)] = [(211)(000)] + ABz * [(210)(000)]
    a3b1[ 6] = a4b0[ 2] + x * a3b0[ 2];    //[(201)(100)] = [(301)(000)] + ABx * [(201)(000)]
    a3b1[ 7] = a4b0[ 4] + y * a3b0[ 2];    //[(201)(010)] = [(211)(000)] + ABy * [(201)(000)]
    a3b1[ 8] = a4b0[ 5] + z * a3b0[ 2];    //[(201)(001)] = [(202)(000)] + ABz * [(201)(000)]
    a3b1[ 9] = a4b0[ 3] + x * a3b0[ 3];    //[(120)(100)] = [(220)(000)] + ABx * [(120)(000)]
    a3b1[10] = a4b0[ 6] + y * a3b0[ 3];    //[(120)(010)] = [(130)(000)] + ABy * [(120)(000)]
    a3b1[11] = a4b0[ 7] + z * a3b0[ 3];    //[(120)(001)] = [(121)(000)] + ABz * [(120)(000)]
    a3b1[12] = a4b0[ 4] + x * a3b0[ 4];    //[(111)(100)] = [(211)(000)] + ABx * [(111)(000)]
    a3b1[13] = a4b0[ 7] + y * a3b0[ 4];    //[(111)(010)] = [(121)(000)] + ABy * [(111)(000)]
    a3b1[14] = a4b0[ 8] + z * a3b0[ 4];    //[(111)(001)] = [(112)(000)] + ABz * [(111)(000)]
    a3b1[15] = a4b0[ 5] + x * a3b0[ 5];    //[(102)(100)] = [(202)(000)] + ABx * [(102)(000)]
    a3b1[16] = a4b0[ 8] + y * a3b0[ 5];    //[(102)(010)] = [(112)(000)] + ABy * [(102)(000)]
    a3b1[17] = a4b0[ 9] + z * a3b0[ 5];    //[(102)(001)] = [(103)(000)] + ABz * [(102)(000)]
    a3b1[18] = a4b0[ 6] + x * a3b0[ 6];    //[(030)(100)] = [(130)(000)] + ABx * [(030)(000)]
    a3b1[19] = a4b0[10] + y * a3b0[ 6];    //[(030)(010)] = [(040)(000)] + ABy * [(030)(000)]
    a3b1[20] = a4b0[11] + z * a3b0[ 6];    //[(030)(001)] = [(031)(000)] + ABz * [(030)(000)]
    a3b1[21] = a4b0[ 7] + x * a3b0[ 7];    //[(021)(100)] = [(121)(000)] + ABx * [(021)(000)]
    a3b1[22] = a4b0[11] + y * a3b0[ 7];    //[(021)(010)] = [(031)(000)] + ABy * [(021)(000)]
    a3b1[23] = a4b0[12] + z * a3b0[ 7];    //[(021)(001)] = [(022)(000)] + ABz * [(021)(000)]
    a3b1[24] = a4b0[ 8] + x * a3b0[ 8];    //[(012)(100)] = [(112)(000)] + ABx * [(012)(000)]
    a3b1[25] = a4b0[12] + y * a3b0[ 8];    //[(012)(010)] = [(022)(000)] + ABy * [(012)(000)]
    a3b1[26] = a4b0[13] + z * a3b0[ 8];    //[(012)(001)] = [(013)(000)] + ABz * [(012)(000)]
    a3b1[27] = a4b0[ 9] + x * a3b0[ 9];    //[(003)(100)] = [(103)(000)] + ABx * [(003)(000)]
    a3b1[28] = a4b0[13] + y * a3b0[ 9];    //[(003)(010)] = [(013)(000)] + ABy * [(003)(000)]
    a3b1[29] = a4b0[14] + z * a3b0[ 9];    //[(003)(001)] = [(004)(000)] + ABz * [(003)(000)]
    a4b1[ 0] = a5b0[ 0] + x * a4b0[ 0];    //[(400)(100)] = [(500)(000)] + ABx * [(400)(000)]
    a4b1[ 3] = a5b0[ 1] + x * a4b0[ 1];    //[(310)(100)] = [(410)(000)] + ABx * [(310)(000)]
    a4b1[ 4] = a5b0[ 3] + y * a4b0[ 1];    //[(310)(010)] = [(320)(000)] + ABy * [(310)(000)]
    a4b1[ 5] = a5b0[ 4] + z * a4b0[ 1];    //[(310)(001)] = [(311)(000)] + ABz * [(310)(000)]
    a4b1[ 6] = a5b0[ 2] + x * a4b0[ 2];    //[(301)(100)] = [(401)(000)] + ABx * [(301)(000)]
    a4b1[ 8] = a5b0[ 5] + z * a4b0[ 2];    //[(301)(001)] = [(302)(000)] + ABz * [(301)(000)]
    a4b1[ 9] = a5b0[ 3] + x * a4b0[ 3];    //[(220)(100)] = [(320)(000)] + ABx * [(220)(000)]
    a4b1[10] = a5b0[ 6] + y * a4b0[ 3];    //[(220)(010)] = [(230)(000)] + ABy * [(220)(000)]
    a4b1[12] = a5b0[ 4] + x * a4b0[ 4];    //[(211)(100)] = [(311)(000)] + ABx * [(211)(000)]
    a4b1[13] = a5b0[ 7] + y * a4b0[ 4];    //[(211)(010)] = [(221)(000)] + ABy * [(211)(000)]
    a4b1[14] = a5b0[ 8] + z * a4b0[ 4];    //[(211)(001)] = [(212)(000)] + ABz * [(211)(000)]
    a4b1[15] = a5b0[ 5] + x * a4b0[ 5];    //[(202)(100)] = [(302)(000)] + ABx * [(202)(000)]
    a4b1[17] = a5b0[ 9] + z * a4b0[ 5];    //[(202)(001)] = [(203)(000)] + ABz * [(202)(000)]
    a4b1[18] = a5b0[ 6] + x * a4b0[ 6];    //[(130)(100)] = [(230)(000)] + ABx * [(130)(000)]
    a4b1[19] = a5b0[10] + y * a4b0[ 6];    //[(130)(010)] = [(140)(000)] + ABy * [(130)(000)]
    a4b1[21] = a5b0[ 7] + x * a4b0[ 7];    //[(121)(100)] = [(221)(000)] + ABx * [(121)(000)]
    a4b1[22] = a5b0[11] + y * a4b0[ 7];    //[(121)(010)] = [(131)(000)] + ABy * [(121)(000)]
    a4b1[23] = a5b0[12] + z * a4b0[ 7];    //[(121)(001)] = [(122)(000)] + ABz * [(121)(000)]
    a4b1[24] = a5b0[ 8] + x * a4b0[ 8];    //[(112)(100)] = [(212)(000)] + ABx * [(112)(000)]
    a4b1[25] = a5b0[12] + y * a4b0[ 8];    //[(112)(010)] = [(122)(000)] + ABy * [(112)(000)]
    a4b1[26] = a5b0[13] + z * a4b0[ 8];    //[(112)(001)] = [(113)(000)] + ABz * [(112)(000)]
    a4b1[27] = a5b0[ 9] + x * a4b0[ 9];    //[(103)(100)] = [(203)(000)] + ABx * [(103)(000)]
    a4b1[28] = a5b0[13] + y * a4b0[ 9];    //[(103)(010)] = [(113)(000)] + ABy * [(103)(000)]
    a4b1[29] = a5b0[14] + z * a4b0[ 9];    //[(103)(001)] = [(104)(000)] + ABz * [(103)(000)]
    a4b1[31] = a5b0[15] + y * a4b0[10];    //[(040)(010)] = [(050)(000)] + ABy * [(040)(000)]
    a4b1[33] = a5b0[11] + x * a4b0[11];    //[(031)(100)] = [(131)(000)] + ABx * [(031)(000)]
    a4b1[34] = a5b0[16] + y * a4b0[11];    //[(031)(010)] = [(041)(000)] + ABy * [(031)(000)]
    a4b1[35] = a5b0[17] + z * a4b0[11];    //[(031)(001)] = [(032)(000)] + ABz * [(031)(000)]
    a4b1[37] = a5b0[17] + y * a4b0[12];    //[(022)(010)] = [(032)(000)] + ABy * [(022)(000)]
    a4b1[38] = a5b0[18] + z * a4b0[12];    //[(022)(001)] = [(023)(000)] + ABz * [(022)(000)]
    a4b1[40] = a5b0[18] + y * a4b0[13];    //[(013)(010)] = [(023)(000)] + ABy * [(013)(000)]
    a4b1[41] = a5b0[19] + z * a4b0[13];    //[(013)(001)] = [(014)(000)] + ABz * [(013)(000)]
    a4b1[44] = a5b0[20] + z * a4b0[14];    //[(004)(001)] = [(005)(000)] + ABz * [(004)(000)]
    a3b2[ 0] = a4b1[ 0] + x * a3b1[ 0];    //[(300)(200)] = [(400)(100)] + ABx * [(300)(100)]
    a3b2[ 1] = a4b1[ 3] + y * a3b1[ 0];    //[(300)(110)] = [(310)(100)] + ABy * [(300)(100)]
    a3b2[ 2] = a4b1[ 6] + z * a3b1[ 0];    //[(300)(101)] = [(301)(100)] + ABz * [(300)(100)]
    a3b2[ 3] = a4b1[ 4] + y * a3b1[ 1];    //[(300)(020)] = [(310)(010)] + ABy * [(300)(010)]
    a3b2[ 4] = a4b1[ 5] + y * a3b1[ 2];    //[(300)(011)] = [(310)(001)] + ABy * [(300)(001)]
    a3b2[ 5] = a4b1[ 8] + z * a3b1[ 2];    //[(300)(002)] = [(301)(001)] + ABz * [(300)(001)]
    a3b2[ 6] = a4b1[ 3] + x * a3b1[ 3];    //[(210)(200)] = [(310)(100)] + ABx * [(210)(100)]
    a3b2[ 7] = a4b1[ 4] + x * a3b1[ 4];    //[(210)(110)] = [(310)(010)] + ABx * [(210)(010)]
    a3b2[ 8] = a4b1[ 5] + x * a3b1[ 5];    //[(210)(101)] = [(310)(001)] + ABx * [(210)(001)]
    a3b2[ 9] = a4b1[10] + y * a3b1[ 4];    //[(210)(020)] = [(220)(010)] + ABy * [(210)(010)]
    a3b2[10] = a4b1[13] + z * a3b1[ 4];    //[(210)(011)] = [(211)(010)] + ABz * [(210)(010)]
    a3b2[11] = a4b1[14] + z * a3b1[ 5];    //[(210)(002)] = [(211)(001)] + ABz * [(210)(001)]
    a3b2[12] = a4b1[ 6] + x * a3b1[ 6];    //[(201)(200)] = [(301)(100)] + ABx * [(201)(100)]
    a3b2[13] = a4b1[12] + y * a3b1[ 6];    //[(201)(110)] = [(211)(100)] + ABy * [(201)(100)]
    a3b2[14] = a4b1[ 8] + x * a3b1[ 8];    //[(201)(101)] = [(301)(001)] + ABx * [(201)(001)]
    a3b2[15] = a4b1[13] + y * a3b1[ 7];    //[(201)(020)] = [(211)(010)] + ABy * [(201)(010)]
    a3b2[16] = a4b1[14] + y * a3b1[ 8];    //[(201)(011)] = [(211)(001)] + ABy * [(201)(001)]
    a3b2[17] = a4b1[17] + z * a3b1[ 8];    //[(201)(002)] = [(202)(001)] + ABz * [(201)(001)]
    a3b2[18] = a4b1[ 9] + x * a3b1[ 9];    //[(120)(200)] = [(220)(100)] + ABx * [(120)(100)]
    a3b2[19] = a4b1[10] + x * a3b1[10];    //[(120)(110)] = [(220)(010)] + ABx * [(120)(010)]
    a3b2[20] = a4b1[21] + z * a3b1[ 9];    //[(120)(101)] = [(121)(100)] + ABz * [(120)(100)]
    a3b2[21] = a4b1[19] + y * a3b1[10];    //[(120)(020)] = [(130)(010)] + ABy * [(120)(010)]
    a3b2[22] = a4b1[22] + z * a3b1[10];    //[(120)(011)] = [(121)(010)] + ABz * [(120)(010)]
    a3b2[23] = a4b1[23] + z * a3b1[11];    //[(120)(002)] = [(121)(001)] + ABz * [(120)(001)]
    a3b2[24] = a4b1[12] + x * a3b1[12];    //[(111)(200)] = [(211)(100)] + ABx * [(111)(100)]
    a3b2[25] = a4b1[13] + x * a3b1[13];    //[(111)(110)] = [(211)(010)] + ABx * [(111)(010)]
    a3b2[26] = a4b1[14] + x * a3b1[14];    //[(111)(101)] = [(211)(001)] + ABx * [(111)(001)]
    a3b2[27] = a4b1[22] + y * a3b1[13];    //[(111)(020)] = [(121)(010)] + ABy * [(111)(010)]
    a3b2[28] = a4b1[23] + y * a3b1[14];    //[(111)(011)] = [(121)(001)] + ABy * [(111)(001)]
    a3b2[29] = a4b1[26] + z * a3b1[14];    //[(111)(002)] = [(112)(001)] + ABz * [(111)(001)]
    a3b2[30] = a4b1[15] + x * a3b1[15];    //[(102)(200)] = [(202)(100)] + ABx * [(102)(100)]
    a3b2[31] = a4b1[24] + y * a3b1[15];    //[(102)(110)] = [(112)(100)] + ABy * [(102)(100)]
    a3b2[32] = a4b1[17] + x * a3b1[17];    //[(102)(101)] = [(202)(001)] + ABx * [(102)(001)]
    a3b2[33] = a4b1[25] + y * a3b1[16];    //[(102)(020)] = [(112)(010)] + ABy * [(102)(010)]
    a3b2[34] = a4b1[26] + y * a3b1[17];    //[(102)(011)] = [(112)(001)] + ABy * [(102)(001)]
    a3b2[35] = a4b1[29] + z * a3b1[17];    //[(102)(002)] = [(103)(001)] + ABz * [(102)(001)]
    a3b2[36] = a4b1[18] + x * a3b1[18];    //[(030)(200)] = [(130)(100)] + ABx * [(030)(100)]
    a3b2[37] = a4b1[19] + x * a3b1[19];    //[(030)(110)] = [(130)(010)] + ABx * [(030)(010)]
    a3b2[38] = a4b1[33] + z * a3b1[18];    //[(030)(101)] = [(031)(100)] + ABz * [(030)(100)]
    a3b2[39] = a4b1[31] + y * a3b1[19];    //[(030)(020)] = [(040)(010)] + ABy * [(030)(010)]
    a3b2[40] = a4b1[34] + z * a3b1[19];    //[(030)(011)] = [(031)(010)] + ABz * [(030)(010)]
    a3b2[41] = a4b1[35] + z * a3b1[20];    //[(030)(002)] = [(031)(001)] + ABz * [(030)(001)]
    a3b2[42] = a4b1[21] + x * a3b1[21];    //[(021)(200)] = [(121)(100)] + ABx * [(021)(100)]
    a3b2[43] = a4b1[22] + x * a3b1[22];    //[(021)(110)] = [(121)(010)] + ABx * [(021)(010)]
    a3b2[44] = a4b1[23] + x * a3b1[23];    //[(021)(101)] = [(121)(001)] + ABx * [(021)(001)]
    a3b2[45] = a4b1[34] + y * a3b1[22];    //[(021)(020)] = [(031)(010)] + ABy * [(021)(010)]
    a3b2[46] = a4b1[35] + y * a3b1[23];    //[(021)(011)] = [(031)(001)] + ABy * [(021)(001)]
    a3b2[47] = a4b1[38] + z * a3b1[23];    //[(021)(002)] = [(022)(001)] + ABz * [(021)(001)]
    a3b2[48] = a4b1[24] + x * a3b1[24];    //[(012)(200)] = [(112)(100)] + ABx * [(012)(100)]
    a3b2[49] = a4b1[25] + x * a3b1[25];    //[(012)(110)] = [(112)(010)] + ABx * [(012)(010)]
    a3b2[50] = a4b1[26] + x * a3b1[26];    //[(012)(101)] = [(112)(001)] + ABx * [(012)(001)]
    a3b2[51] = a4b1[37] + y * a3b1[25];    //[(012)(020)] = [(022)(010)] + ABy * [(012)(010)]
    a3b2[52] = a4b1[38] + y * a3b1[26];    //[(012)(011)] = [(022)(001)] + ABy * [(012)(001)]
    a3b2[53] = a4b1[41] + z * a3b1[26];    //[(012)(002)] = [(013)(001)] + ABz * [(012)(001)]
    a3b2[54] = a4b1[27] + x * a3b1[27];    //[(003)(200)] = [(103)(100)] + ABx * [(003)(100)]
    a3b2[55] = a4b1[28] + x * a3b1[28];    //[(003)(110)] = [(103)(010)] + ABx * [(003)(010)]
    a3b2[56] = a4b1[29] + x * a3b1[29];    //[(003)(101)] = [(103)(001)] + ABx * [(003)(001)]
    a3b2[57] = a4b1[40] + y * a3b1[28];    //[(003)(020)] = [(013)(010)] + ABy * [(003)(010)]
    a3b2[58] = a4b1[41] + y * a3b1[29];    //[(003)(011)] = [(013)(001)] + ABy * [(003)(001)]
    a3b2[59] = a4b1[44] + z * a3b1[29];    //[(003)(002)] = [(004)(001)] + ABz * [(003)(001)]
}

} // namespace (hgp)
} // namespace (nhfInt)
