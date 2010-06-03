
        // 
        // This file was automatically generated by a sage script,
        // which should be in the file w_coefficient.sage, and it
        // probably should not be edited manually. Instead, edit
        // that file, and rerun the script.
        //
        // To run the script, you can, for example, run sage from
        // the command line, attach the file w_coefficient.sage,
        // and then run the function write_w_coefficient_file(limit)
        // with limit equal to the one more than the largest j you
        // wish these coefficients to be precomputed with.
        //
    const int A_range = 11;
Complex A[11][11] = {{Complex(0.50000000000000000,0.50000000000000000),
Complex(-0.14104739588693907,0.00000000000000000),
Complex(0.039788735772973836,-0.039788735772973836),
Complex(0.00000000000000000,0.033672585398468728),
Complex(-0.018997721932938333,-0.018997721932938333),
Complex(0.026795792064251384,0.00000000000000000),
Complex(-0.022676860148343390,0.022676860148343390),
Complex(0.00000000000000000,-0.044779168991425991),
Complex(0.050527881409774461,0.050527881409774461),
Complex(-0.12828286966558983,0.00000000000000000),
Complex(0.18093964703235055,-0.18093964703235055)},
{Complex(0.0,0.0),
Complex(0.25000000000000000,0.25000000000000000),
Complex(-0.14104739588693907,0.00000000000000000),
Complex(0.059683103659460751,-0.059683103659460751),
Complex(0.00000000000000000,0.067345170796937456),
Complex(-0.047494304832345832,-0.047494304832345832),
Complex(0.080387376192754142,0.00000000000000000),
Complex(-0.079369010519201866,0.079369010519201866),
Complex(0.00000000000000000,-0.17911667596570396),
Complex(0.22737546634398509,0.22737546634398509),
Complex(-0.64141434832794908,0.00000000000000000)},
{Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.12500000000000000,0.12500000000000000),
Complex(-0.10578554691520431,0.00000000000000000),
Complex(0.059683103659460751,-0.059683103659460751),
Complex(0.00000000000000000,0.084181463496171824),
Complex(-0.071241457248518741,-0.071241457248518741),
Complex(0.14067790833731977,0.00000000000000000),
Complex(-0.15873802103840373,0.15873802103840373),
Complex(0.00000000000000000,-0.40301252092283391),
Complex(0.56843866585996272,0.56843866585996272)},
{Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.062500000000000000,0.062500000000000000),
Complex(-0.070523697943469535,0.00000000000000000),
Complex(0.049735919716217290,-0.049735919716217290),
Complex(0.00000000000000000,0.084181463496171824),
Complex(-0.083115033456605203,-0.083115033456605203),
Complex(0.18757054444975968,0.00000000000000000),
Complex(-0.23810703155760560,0.23810703155760560),
Complex(0.00000000000000000,-0.67168753487138988)},
{Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.031250000000000000,0.031250000000000000),
Complex(-0.044077311214668458,0.00000000000000000),
Complex(0.037301939787162966,-0.037301939787162966),
Complex(0.00000000000000000,0.073658780559150344),
Complex(-0.083115033456605203,-0.083115033456605203),
Complex(0.21101686250597965,0.00000000000000000),
Complex(-0.29763378944700702,0.29763378944700702)},
{Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.015625000000000000,0.015625000000000000),
Complex(-0.026446386728801077,0.00000000000000000),
Complex(0.026111357851014077,-0.026111357851014077),
Complex(0.00000000000000000,0.058927024447320279),
Complex(-0.074803530110944677,-0.074803530110944677),
Complex(0.21101686250597965,0.00000000000000000)},
{Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0078125000000000000,0.0078125000000000000),
Complex(-0.015427058925133961,0.00000000000000000),
Complex(0.017407571900676051,-0.017407571900676051),
Complex(0.00000000000000000,0.044195268335490208),
Complex(-0.062336275092453902,-0.062336275092453902)},
{Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0039062500000000000,0.0039062500000000000),
Complex(-0.0088154622429336919,0.00000000000000000),
Complex(0.011190581936148891,-0.011190581936148891),
Complex(0.00000000000000000,0.031568048811064432)},
{Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0019531250000000000,0.0019531250000000000),
Complex(-0.0049586975116502020,0.00000000000000000),
Complex(0.0069941137100930570,-0.0069941137100930570)},
{Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.00097656250000000000,0.00097656250000000000),
Complex(-0.0027548319509167786,0.00000000000000000)},
{Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.00048828125000000000,0.00048828125000000000)}};
const int B_range = 11;
Complex B[11][11][11] = {{{Complex(1.0000000000000000,0.00000000000000000),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0)},
{Complex(0.0,0.0),
Complex(1.7724538509055161,1.7724538509055161),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0)},
{Complex(-1.0000000000000000,0.00000000000000000),
Complex(0.0,0.0),
Complex(0.00000000000000000,3.1415926535897931),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0)},
{Complex(0.0,0.0),
Complex(-1.7724538509055161,-1.7724538509055161),
Complex(0.0,0.0),
Complex(-1.8561093322772360,1.8561093322772360),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0)},
{Complex(0.50000000000000000,0.00000000000000000),
Complex(0.0,0.0),
Complex(0.00000000000000000,-3.1415926535897931),
Complex(0.0,0.0),
Complex(-1.6449340668482264,0.00000000000000000),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0)},
{Complex(0.0,0.0),
Complex(0.88622692545275805,0.88622692545275805),
Complex(0.0,0.0),
Complex(1.8561093322772360,-1.8561093322772360),
Complex(0.0,0.0),
Complex(-0.58311394425416208,-0.58311394425416208),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0)},
{Complex(-0.16666666666666666,0.00000000000000000),
Complex(0.0,0.0),
Complex(0.00000000000000000,1.5707963267948966),
Complex(0.0,0.0),
Complex(1.6449340668482264,0.00000000000000000),
Complex(0.0,0.0),
Complex(0.00000000000000000,-0.34451418533666467),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0)},
{Complex(0.0,0.0),
Complex(-0.29540897515091935,-0.29540897515091935),
Complex(0.0,0.0),
Complex(-0.92805466613861798,0.92805466613861798),
Complex(0.0,0.0),
Complex(0.58311394425416208,0.58311394425416208),
Complex(0.0,0.0),
Complex(0.087233642070221135,-0.087233642070221135),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0)},
{Complex(0.041666666666666664,0.00000000000000000),
Complex(0.0,0.0),
Complex(0.00000000000000000,-0.52359877559829893),
Complex(0.0,0.0),
Complex(-0.82246703342411320,0.00000000000000000),
Complex(0.0,0.0),
Complex(0.00000000000000000,0.34451418533666467),
Complex(0.0,0.0),
Complex(0.038654401203969221,0.00000000000000000),
Complex(0.0,0.0),
Complex(0.0,0.0)},
{Complex(0.0,0.0),
Complex(0.073852243787729838,0.073852243787729838),
Complex(0.0,0.0),
Complex(0.30935155537953934,-0.30935155537953934),
Complex(0.0,0.0),
Complex(-0.29155697212708104,-0.29155697212708104),
Complex(0.0,0.0),
Complex(-0.087233642070221135,0.087233642070221135),
Complex(0.0,0.0),
Complex(0.0076125713631580065,0.0076125713631580065),
Complex(0.0,0.0)},
{Complex(-0.0083333333333333332,0.00000000000000000),
Complex(0.0,0.0),
Complex(0.00000000000000000,0.13089969389957473),
Complex(0.0,0.0),
Complex(0.27415567780803773,0.00000000000000000),
Complex(0.0,0.0),
Complex(0.00000000000000000,-0.17225709266833233),
Complex(0.0,0.0),
Complex(-0.038654401203969221,0.00000000000000000),
Complex(0.0,0.0),
Complex(0.00000000000000000,0.0026985862855844925)}},
{{Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0)},
{Complex(1.0000000000000000,0.00000000000000000),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0)},
{Complex(0.0,0.0),
Complex(1.7724538509055161,1.7724538509055161),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0)},
{Complex(-1.0000000000000000,0.00000000000000000),
Complex(0.0,0.0),
Complex(0.00000000000000000,3.1415926535897931),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0)},
{Complex(0.0,0.0),
Complex(-1.7724538509055161,-1.7724538509055161),
Complex(0.0,0.0),
Complex(-1.8561093322772360,1.8561093322772360),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0)},
{Complex(0.50000000000000000,0.00000000000000000),
Complex(0.0,0.0),
Complex(0.00000000000000000,-3.1415926535897931),
Complex(0.0,0.0),
Complex(-1.6449340668482264,0.00000000000000000),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0)},
{Complex(0.0,0.0),
Complex(0.88622692545275805,0.88622692545275805),
Complex(0.0,0.0),
Complex(1.8561093322772360,-1.8561093322772360),
Complex(0.0,0.0),
Complex(-0.58311394425416208,-0.58311394425416208),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0)},
{Complex(-0.16666666666666666,0.00000000000000000),
Complex(0.0,0.0),
Complex(0.00000000000000000,1.5707963267948966),
Complex(0.0,0.0),
Complex(1.6449340668482264,0.00000000000000000),
Complex(0.0,0.0),
Complex(0.00000000000000000,-0.34451418533666467),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0)},
{Complex(0.0,0.0),
Complex(-0.29540897515091935,-0.29540897515091935),
Complex(0.0,0.0),
Complex(-0.92805466613861798,0.92805466613861798),
Complex(0.0,0.0),
Complex(0.58311394425416208,0.58311394425416208),
Complex(0.0,0.0),
Complex(0.087233642070221135,-0.087233642070221135),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0)},
{Complex(0.041666666666666664,0.00000000000000000),
Complex(0.0,0.0),
Complex(0.00000000000000000,-0.52359877559829893),
Complex(0.0,0.0),
Complex(-0.82246703342411320,0.00000000000000000),
Complex(0.0,0.0),
Complex(0.00000000000000000,0.34451418533666467),
Complex(0.0,0.0),
Complex(0.038654401203969221,0.00000000000000000),
Complex(0.0,0.0),
Complex(0.0,0.0)},
{Complex(0.0,0.0),
Complex(0.073852243787729838,0.073852243787729838),
Complex(0.0,0.0),
Complex(0.30935155537953934,-0.30935155537953934),
Complex(0.0,0.0),
Complex(-0.29155697212708104,-0.29155697212708104),
Complex(0.0,0.0),
Complex(-0.087233642070221135,0.087233642070221135),
Complex(0.0,0.0),
Complex(0.0076125713631580065,0.0076125713631580065),
Complex(0.0,0.0)}},
{{Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0)},
{Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0)},
{Complex(1.0000000000000000,0.00000000000000000),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0)},
{Complex(0.0,0.0),
Complex(1.7724538509055161,1.7724538509055161),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0)},
{Complex(-1.0000000000000000,0.00000000000000000),
Complex(0.0,0.0),
Complex(0.00000000000000000,3.1415926535897931),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0)},
{Complex(0.0,0.0),
Complex(-1.7724538509055161,-1.7724538509055161),
Complex(0.0,0.0),
Complex(-1.8561093322772360,1.8561093322772360),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0)},
{Complex(0.50000000000000000,0.00000000000000000),
Complex(0.0,0.0),
Complex(0.00000000000000000,-3.1415926535897931),
Complex(0.0,0.0),
Complex(-1.6449340668482264,0.00000000000000000),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0)},
{Complex(0.0,0.0),
Complex(0.88622692545275805,0.88622692545275805),
Complex(0.0,0.0),
Complex(1.8561093322772360,-1.8561093322772360),
Complex(0.0,0.0),
Complex(-0.58311394425416208,-0.58311394425416208),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0)},
{Complex(-0.16666666666666666,0.00000000000000000),
Complex(0.0,0.0),
Complex(0.00000000000000000,1.5707963267948966),
Complex(0.0,0.0),
Complex(1.6449340668482264,0.00000000000000000),
Complex(0.0,0.0),
Complex(0.00000000000000000,-0.34451418533666467),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0)},
{Complex(0.0,0.0),
Complex(-0.29540897515091935,-0.29540897515091935),
Complex(0.0,0.0),
Complex(-0.92805466613861798,0.92805466613861798),
Complex(0.0,0.0),
Complex(0.58311394425416208,0.58311394425416208),
Complex(0.0,0.0),
Complex(0.087233642070221135,-0.087233642070221135),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0)},
{Complex(0.041666666666666664,0.00000000000000000),
Complex(0.0,0.0),
Complex(0.00000000000000000,-0.52359877559829893),
Complex(0.0,0.0),
Complex(-0.82246703342411320,0.00000000000000000),
Complex(0.0,0.0),
Complex(0.00000000000000000,0.34451418533666467),
Complex(0.0,0.0),
Complex(0.038654401203969221,0.00000000000000000),
Complex(0.0,0.0),
Complex(0.0,0.0)}},
{{Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0)},
{Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0)},
{Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0)},
{Complex(1.0000000000000000,0.00000000000000000),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0)},
{Complex(0.0,0.0),
Complex(1.7724538509055161,1.7724538509055161),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0)},
{Complex(-1.0000000000000000,0.00000000000000000),
Complex(0.0,0.0),
Complex(0.00000000000000000,3.1415926535897931),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0)},
{Complex(0.0,0.0),
Complex(-1.7724538509055161,-1.7724538509055161),
Complex(0.0,0.0),
Complex(-1.8561093322772360,1.8561093322772360),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0)},
{Complex(0.50000000000000000,0.00000000000000000),
Complex(0.0,0.0),
Complex(0.00000000000000000,-3.1415926535897931),
Complex(0.0,0.0),
Complex(-1.6449340668482264,0.00000000000000000),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0)},
{Complex(0.0,0.0),
Complex(0.88622692545275805,0.88622692545275805),
Complex(0.0,0.0),
Complex(1.8561093322772360,-1.8561093322772360),
Complex(0.0,0.0),
Complex(-0.58311394425416208,-0.58311394425416208),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0)},
{Complex(-0.16666666666666666,0.00000000000000000),
Complex(0.0,0.0),
Complex(0.00000000000000000,1.5707963267948966),
Complex(0.0,0.0),
Complex(1.6449340668482264,0.00000000000000000),
Complex(0.0,0.0),
Complex(0.00000000000000000,-0.34451418533666467),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0)},
{Complex(0.0,0.0),
Complex(-0.29540897515091935,-0.29540897515091935),
Complex(0.0,0.0),
Complex(-0.92805466613861798,0.92805466613861798),
Complex(0.0,0.0),
Complex(0.58311394425416208,0.58311394425416208),
Complex(0.0,0.0),
Complex(0.087233642070221135,-0.087233642070221135),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0)}},
{{Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0)},
{Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0)},
{Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0)},
{Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0)},
{Complex(1.0000000000000000,0.00000000000000000),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0)},
{Complex(0.0,0.0),
Complex(1.7724538509055161,1.7724538509055161),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0)},
{Complex(-1.0000000000000000,0.00000000000000000),
Complex(0.0,0.0),
Complex(0.00000000000000000,3.1415926535897931),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0)},
{Complex(0.0,0.0),
Complex(-1.7724538509055161,-1.7724538509055161),
Complex(0.0,0.0),
Complex(-1.8561093322772360,1.8561093322772360),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0)},
{Complex(0.50000000000000000,0.00000000000000000),
Complex(0.0,0.0),
Complex(0.00000000000000000,-3.1415926535897931),
Complex(0.0,0.0),
Complex(-1.6449340668482264,0.00000000000000000),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0)},
{Complex(0.0,0.0),
Complex(0.88622692545275805,0.88622692545275805),
Complex(0.0,0.0),
Complex(1.8561093322772360,-1.8561093322772360),
Complex(0.0,0.0),
Complex(-0.58311394425416208,-0.58311394425416208),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0)},
{Complex(-0.16666666666666666,0.00000000000000000),
Complex(0.0,0.0),
Complex(0.00000000000000000,1.5707963267948966),
Complex(0.0,0.0),
Complex(1.6449340668482264,0.00000000000000000),
Complex(0.0,0.0),
Complex(0.00000000000000000,-0.34451418533666467),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0)}},
{{Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0)},
{Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0)},
{Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0)},
{Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0)},
{Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0)},
{Complex(1.0000000000000000,0.00000000000000000),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0)},
{Complex(0.0,0.0),
Complex(1.7724538509055161,1.7724538509055161),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0)},
{Complex(-1.0000000000000000,0.00000000000000000),
Complex(0.0,0.0),
Complex(0.00000000000000000,3.1415926535897931),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0)},
{Complex(0.0,0.0),
Complex(-1.7724538509055161,-1.7724538509055161),
Complex(0.0,0.0),
Complex(-1.8561093322772360,1.8561093322772360),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0)},
{Complex(0.50000000000000000,0.00000000000000000),
Complex(0.0,0.0),
Complex(0.00000000000000000,-3.1415926535897931),
Complex(0.0,0.0),
Complex(-1.6449340668482264,0.00000000000000000),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0)},
{Complex(0.0,0.0),
Complex(0.88622692545275805,0.88622692545275805),
Complex(0.0,0.0),
Complex(1.8561093322772360,-1.8561093322772360),
Complex(0.0,0.0),
Complex(-0.58311394425416208,-0.58311394425416208),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0)}},
{{Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0)},
{Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0)},
{Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0)},
{Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0)},
{Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0)},
{Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0)},
{Complex(1.0000000000000000,0.00000000000000000),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0)},
{Complex(0.0,0.0),
Complex(1.7724538509055161,1.7724538509055161),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0)},
{Complex(-1.0000000000000000,0.00000000000000000),
Complex(0.0,0.0),
Complex(0.00000000000000000,3.1415926535897931),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0)},
{Complex(0.0,0.0),
Complex(-1.7724538509055161,-1.7724538509055161),
Complex(0.0,0.0),
Complex(-1.8561093322772360,1.8561093322772360),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0)},
{Complex(0.50000000000000000,0.00000000000000000),
Complex(0.0,0.0),
Complex(0.00000000000000000,-3.1415926535897931),
Complex(0.0,0.0),
Complex(-1.6449340668482264,0.00000000000000000),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0)}},
{{Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0)},
{Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0)},
{Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0)},
{Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0)},
{Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0)},
{Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0)},
{Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0)},
{Complex(1.0000000000000000,0.00000000000000000),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0)},
{Complex(0.0,0.0),
Complex(1.7724538509055161,1.7724538509055161),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0)},
{Complex(-1.0000000000000000,0.00000000000000000),
Complex(0.0,0.0),
Complex(0.00000000000000000,3.1415926535897931),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0)},
{Complex(0.0,0.0),
Complex(-1.7724538509055161,-1.7724538509055161),
Complex(0.0,0.0),
Complex(-1.8561093322772360,1.8561093322772360),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0)}},
{{Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0)},
{Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0)},
{Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0)},
{Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0)},
{Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0)},
{Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0)},
{Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0)},
{Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0)},
{Complex(1.0000000000000000,0.00000000000000000),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0)},
{Complex(0.0,0.0),
Complex(1.7724538509055161,1.7724538509055161),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0)},
{Complex(-1.0000000000000000,0.00000000000000000),
Complex(0.0,0.0),
Complex(0.00000000000000000,3.1415926535897931),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0)}},
{{Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0)},
{Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0)},
{Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0)},
{Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0)},
{Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0)},
{Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0)},
{Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0)},
{Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0)},
{Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0)},
{Complex(1.0000000000000000,0.00000000000000000),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0)},
{Complex(0.0,0.0),
Complex(1.7724538509055161,1.7724538509055161),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0)}},
{{Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0)},
{Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0)},
{Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0)},
{Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0)},
{Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0)},
{Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0)},
{Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0)},
{Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0)},
{Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0)},
{Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0)},
{Complex(1.0000000000000000,0.00000000000000000),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0),
Complex(0.0,0.0)}}};
