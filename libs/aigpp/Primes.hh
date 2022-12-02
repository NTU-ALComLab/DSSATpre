/**************************************************************
 *
 *       AIGPP Package // Primes.hh
 *
 *	Author:
 *         Florian Pigorsch
 *	  University of Freiburg
 *         pigorsch@informatik.uni-freiburg.de
 *
 *	Last revision:
 *         $Revision: 717 $
 *	  $Author: pigorsch $
 *         $Date$
 *
 ***************************************************************/

#ifndef AIGPP_PRIMES_HH
#define AIGPP_PRIMES_HH

namespace aigpp {
enum
{
    NumPowerPrimes = 28
};

static const unsigned long PowerPrimes[NumPowerPrimes]
    = {53ul,        97ul,        193ul,       389ul,       769ul,        1543ul,       3079ul,
       6151ul,      12289ul,     24593ul,     49157ul,     98317ul,      196613ul,     393241ul,
       786433ul,    1572869ul,   3145739ul,   6291469ul,   12582917ul,   25165843ul,   50331653ul,
       100663319ul, 201326611ul, 402653189ul, 805306457ul, 1610612741ul, 3221225473ul, 4294967291ul};

enum
{
    NumPrimes = 1000
};

static const unsigned long Primes[NumPrimes]
    = {2ul,    3ul,    5ul,    7ul,    11ul,   13ul,   17ul,   19ul,   23ul,   29ul,   31ul,   37ul,   41ul,   43ul,
       47ul,   53ul,   59ul,   61ul,   67ul,   71ul,   73ul,   79ul,   83ul,   89ul,   97ul,   101ul,  103ul,  107ul,
       109ul,  113ul,  127ul,  131ul,  137ul,  139ul,  149ul,  151ul,  157ul,  163ul,  167ul,  173ul,  179ul,  181ul,
       191ul,  193ul,  197ul,  199ul,  211ul,  223ul,  227ul,  229ul,  233ul,  239ul,  241ul,  251ul,  257ul,  263ul,
       269ul,  271ul,  277ul,  281ul,  283ul,  293ul,  307ul,  311ul,  313ul,  317ul,  331ul,  337ul,  347ul,  349ul,
       353ul,  359ul,  367ul,  373ul,  379ul,  383ul,  389ul,  397ul,  401ul,  409ul,  419ul,  421ul,  431ul,  433ul,
       439ul,  443ul,  449ul,  457ul,  461ul,  463ul,  467ul,  479ul,  487ul,  491ul,  499ul,  503ul,  509ul,  521ul,
       523ul,  541ul,  547ul,  557ul,  563ul,  569ul,  571ul,  577ul,  587ul,  593ul,  599ul,  601ul,  607ul,  613ul,
       617ul,  619ul,  631ul,  641ul,  643ul,  647ul,  653ul,  659ul,  661ul,  673ul,  677ul,  683ul,  691ul,  701ul,
       709ul,  719ul,  727ul,  733ul,  739ul,  743ul,  751ul,  757ul,  761ul,  769ul,  773ul,  787ul,  797ul,  809ul,
       811ul,  821ul,  823ul,  827ul,  829ul,  839ul,  853ul,  857ul,  859ul,  863ul,  877ul,  881ul,  883ul,  887ul,
       907ul,  911ul,  919ul,  929ul,  937ul,  941ul,  947ul,  953ul,  967ul,  971ul,  977ul,  983ul,  991ul,  997ul,
       1009ul, 1013ul, 1019ul, 1021ul, 1031ul, 1033ul, 1039ul, 1049ul, 1051ul, 1061ul, 1063ul, 1069ul, 1087ul, 1091ul,
       1093ul, 1097ul, 1103ul, 1109ul, 1117ul, 1123ul, 1129ul, 1151ul, 1153ul, 1163ul, 1171ul, 1181ul, 1187ul, 1193ul,
       1201ul, 1213ul, 1217ul, 1223ul, 1229ul, 1231ul, 1237ul, 1249ul, 1259ul, 1277ul, 1279ul, 1283ul, 1289ul, 1291ul,
       1297ul, 1301ul, 1303ul, 1307ul, 1319ul, 1321ul, 1327ul, 1361ul, 1367ul, 1373ul, 1381ul, 1399ul, 1409ul, 1423ul,
       1427ul, 1429ul, 1433ul, 1439ul, 1447ul, 1451ul, 1453ul, 1459ul, 1471ul, 1481ul, 1483ul, 1487ul, 1489ul, 1493ul,
       1499ul, 1511ul, 1523ul, 1531ul, 1543ul, 1549ul, 1553ul, 1559ul, 1567ul, 1571ul, 1579ul, 1583ul, 1597ul, 1601ul,
       1607ul, 1609ul, 1613ul, 1619ul, 1621ul, 1627ul, 1637ul, 1657ul, 1663ul, 1667ul, 1669ul, 1693ul, 1697ul, 1699ul,
       1709ul, 1721ul, 1723ul, 1733ul, 1741ul, 1747ul, 1753ul, 1759ul, 1777ul, 1783ul, 1787ul, 1789ul, 1801ul, 1811ul,
       1823ul, 1831ul, 1847ul, 1861ul, 1867ul, 1871ul, 1873ul, 1877ul, 1879ul, 1889ul, 1901ul, 1907ul, 1913ul, 1931ul,
       1933ul, 1949ul, 1951ul, 1973ul, 1979ul, 1987ul, 1993ul, 1997ul, 1999ul, 2003ul, 2011ul, 2017ul, 2027ul, 2029ul,
       2039ul, 2053ul, 2063ul, 2069ul, 2081ul, 2083ul, 2087ul, 2089ul, 2099ul, 2111ul, 2113ul, 2129ul, 2131ul, 2137ul,
       2141ul, 2143ul, 2153ul, 2161ul, 2179ul, 2203ul, 2207ul, 2213ul, 2221ul, 2237ul, 2239ul, 2243ul, 2251ul, 2267ul,
       2269ul, 2273ul, 2281ul, 2287ul, 2293ul, 2297ul, 2309ul, 2311ul, 2333ul, 2339ul, 2341ul, 2347ul, 2351ul, 2357ul,
       2371ul, 2377ul, 2381ul, 2383ul, 2389ul, 2393ul, 2399ul, 2411ul, 2417ul, 2423ul, 2437ul, 2441ul, 2447ul, 2459ul,
       2467ul, 2473ul, 2477ul, 2503ul, 2521ul, 2531ul, 2539ul, 2543ul, 2549ul, 2551ul, 2557ul, 2579ul, 2591ul, 2593ul,
       2609ul, 2617ul, 2621ul, 2633ul, 2647ul, 2657ul, 2659ul, 2663ul, 2671ul, 2677ul, 2683ul, 2687ul, 2689ul, 2693ul,
       2699ul, 2707ul, 2711ul, 2713ul, 2719ul, 2729ul, 2731ul, 2741ul, 2749ul, 2753ul, 2767ul, 2777ul, 2789ul, 2791ul,
       2797ul, 2801ul, 2803ul, 2819ul, 2833ul, 2837ul, 2843ul, 2851ul, 2857ul, 2861ul, 2879ul, 2887ul, 2897ul, 2903ul,
       2909ul, 2917ul, 2927ul, 2939ul, 2953ul, 2957ul, 2963ul, 2969ul, 2971ul, 2999ul, 3001ul, 3011ul, 3019ul, 3023ul,
       3037ul, 3041ul, 3049ul, 3061ul, 3067ul, 3079ul, 3083ul, 3089ul, 3109ul, 3119ul, 3121ul, 3137ul, 3163ul, 3167ul,
       3169ul, 3181ul, 3187ul, 3191ul, 3203ul, 3209ul, 3217ul, 3221ul, 3229ul, 3251ul, 3253ul, 3257ul, 3259ul, 3271ul,
       3299ul, 3301ul, 3307ul, 3313ul, 3319ul, 3323ul, 3329ul, 3331ul, 3343ul, 3347ul, 3359ul, 3361ul, 3371ul, 3373ul,
       3389ul, 3391ul, 3407ul, 3413ul, 3433ul, 3449ul, 3457ul, 3461ul, 3463ul, 3467ul, 3469ul, 3491ul, 3499ul, 3511ul,
       3517ul, 3527ul, 3529ul, 3533ul, 3539ul, 3541ul, 3547ul, 3557ul, 3559ul, 3571ul, 3581ul, 3583ul, 3593ul, 3607ul,
       3613ul, 3617ul, 3623ul, 3631ul, 3637ul, 3643ul, 3659ul, 3671ul, 3673ul, 3677ul, 3691ul, 3697ul, 3701ul, 3709ul,
       3719ul, 3727ul, 3733ul, 3739ul, 3761ul, 3767ul, 3769ul, 3779ul, 3793ul, 3797ul, 3803ul, 3821ul, 3823ul, 3833ul,
       3847ul, 3851ul, 3853ul, 3863ul, 3877ul, 3881ul, 3889ul, 3907ul, 3911ul, 3917ul, 3919ul, 3923ul, 3929ul, 3931ul,
       3943ul, 3947ul, 3967ul, 3989ul, 4001ul, 4003ul, 4007ul, 4013ul, 4019ul, 4021ul, 4027ul, 4049ul, 4051ul, 4057ul,
       4073ul, 4079ul, 4091ul, 4093ul, 4099ul, 4111ul, 4127ul, 4129ul, 4133ul, 4139ul, 4153ul, 4157ul, 4159ul, 4177ul,
       4201ul, 4211ul, 4217ul, 4219ul, 4229ul, 4231ul, 4241ul, 4243ul, 4253ul, 4259ul, 4261ul, 4271ul, 4273ul, 4283ul,
       4289ul, 4297ul, 4327ul, 4337ul, 4339ul, 4349ul, 4357ul, 4363ul, 4373ul, 4391ul, 4397ul, 4409ul, 4421ul, 4423ul,
       4441ul, 4447ul, 4451ul, 4457ul, 4463ul, 4481ul, 4483ul, 4493ul, 4507ul, 4513ul, 4517ul, 4519ul, 4523ul, 4547ul,
       4549ul, 4561ul, 4567ul, 4583ul, 4591ul, 4597ul, 4603ul, 4621ul, 4637ul, 4639ul, 4643ul, 4649ul, 4651ul, 4657ul,
       4663ul, 4673ul, 4679ul, 4691ul, 4703ul, 4721ul, 4723ul, 4729ul, 4733ul, 4751ul, 4759ul, 4783ul, 4787ul, 4789ul,
       4793ul, 4799ul, 4801ul, 4813ul, 4817ul, 4831ul, 4861ul, 4871ul, 4877ul, 4889ul, 4903ul, 4909ul, 4919ul, 4931ul,
       4933ul, 4937ul, 4943ul, 4951ul, 4957ul, 4967ul, 4969ul, 4973ul, 4987ul, 4993ul, 4999ul, 5003ul, 5009ul, 5011ul,
       5021ul, 5023ul, 5039ul, 5051ul, 5059ul, 5077ul, 5081ul, 5087ul, 5099ul, 5101ul, 5107ul, 5113ul, 5119ul, 5147ul,
       5153ul, 5167ul, 5171ul, 5179ul, 5189ul, 5197ul, 5209ul, 5227ul, 5231ul, 5233ul, 5237ul, 5261ul, 5273ul, 5279ul,
       5281ul, 5297ul, 5303ul, 5309ul, 5323ul, 5333ul, 5347ul, 5351ul, 5381ul, 5387ul, 5393ul, 5399ul, 5407ul, 5413ul,
       5417ul, 5419ul, 5431ul, 5437ul, 5441ul, 5443ul, 5449ul, 5471ul, 5477ul, 5479ul, 5483ul, 5501ul, 5503ul, 5507ul,
       5519ul, 5521ul, 5527ul, 5531ul, 5557ul, 5563ul, 5569ul, 5573ul, 5581ul, 5591ul, 5623ul, 5639ul, 5641ul, 5647ul,
       5651ul, 5653ul, 5657ul, 5659ul, 5669ul, 5683ul, 5689ul, 5693ul, 5701ul, 5711ul, 5717ul, 5737ul, 5741ul, 5743ul,
       5749ul, 5779ul, 5783ul, 5791ul, 5801ul, 5807ul, 5813ul, 5821ul, 5827ul, 5839ul, 5843ul, 5849ul, 5851ul, 5857ul,
       5861ul, 5867ul, 5869ul, 5879ul, 5881ul, 5897ul, 5903ul, 5923ul, 5927ul, 5939ul, 5953ul, 5981ul, 5987ul, 6007ul,
       6011ul, 6029ul, 6037ul, 6043ul, 6047ul, 6053ul, 6067ul, 6073ul, 6079ul, 6089ul, 6091ul, 6101ul, 6113ul, 6121ul,
       6131ul, 6133ul, 6143ul, 6151ul, 6163ul, 6173ul, 6197ul, 6199ul, 6203ul, 6211ul, 6217ul, 6221ul, 6229ul, 6247ul,
       6257ul, 6263ul, 6269ul, 6271ul, 6277ul, 6287ul, 6299ul, 6301ul, 6311ul, 6317ul, 6323ul, 6329ul, 6337ul, 6343ul,
       6353ul, 6359ul, 6361ul, 6367ul, 6373ul, 6379ul, 6389ul, 6397ul, 6421ul, 6427ul, 6449ul, 6451ul, 6469ul, 6473ul,
       6481ul, 6491ul, 6521ul, 6529ul, 6547ul, 6551ul, 6553ul, 6563ul, 6569ul, 6571ul, 6577ul, 6581ul, 6599ul, 6607ul,
       6619ul, 6637ul, 6653ul, 6659ul, 6661ul, 6673ul, 6679ul, 6689ul, 6691ul, 6701ul, 6703ul, 6709ul, 6719ul, 6733ul,
       6737ul, 6761ul, 6763ul, 6779ul, 6781ul, 6791ul, 6793ul, 6803ul, 6823ul, 6827ul, 6829ul, 6833ul, 6841ul, 6857ul,
       6863ul, 6869ul, 6871ul, 6883ul, 6899ul, 6907ul, 6911ul, 6917ul, 6947ul, 6949ul, 6959ul, 6961ul, 6967ul, 6971ul,
       6977ul, 6983ul, 6991ul, 6997ul, 7001ul, 7013ul, 7019ul, 7027ul, 7039ul, 7043ul, 7057ul, 7069ul, 7079ul, 7103ul,
       7109ul, 7121ul, 7127ul, 7129ul, 7151ul, 7159ul, 7177ul, 7187ul, 7193ul, 7207ul, 7211ul, 7213ul, 7219ul, 7229ul,
       7237ul, 7243ul, 7247ul, 7253ul, 7283ul, 7297ul, 7307ul, 7309ul, 7321ul, 7331ul, 7333ul, 7349ul, 7351ul, 7369ul,
       7393ul, 7411ul, 7417ul, 7433ul, 7451ul, 7457ul, 7459ul, 7477ul, 7481ul, 7487ul, 7489ul, 7499ul, 7507ul, 7517ul,
       7523ul, 7529ul, 7537ul, 7541ul, 7547ul, 7549ul, 7559ul, 7561ul, 7573ul, 7577ul, 7583ul, 7589ul, 7591ul, 7603ul,
       7607ul, 7621ul, 7639ul, 7643ul, 7649ul, 7669ul, 7673ul, 7681ul, 7687ul, 7691ul, 7699ul, 7703ul, 7717ul, 7723ul,
       7727ul, 7741ul, 7753ul, 7757ul, 7759ul, 7789ul, 7793ul, 7817ul, 7823ul, 7829ul, 7841ul, 7853ul, 7867ul, 7873ul,
       7877ul, 7879ul, 7883ul, 7901ul, 7907ul, 7919ul};

}  // namespace aigpp

#endif /* AIGPP_PRIMES_HH */
