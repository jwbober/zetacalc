#include <queue>
#include <fstream>

#include "theta_sums.h"
#include "main_sum.h"

#include "log.h"

using namespace std;

Complex zeta_block_mpfr(mpfr_t v, unsigned int K, mpfr_t t) {
    mpz_t vv;
    mpz_init(vv);

    mpfr_get_z(vv, v, GMP_RNDN);

    Complex S = zeta_block_mpfr(vv, K, t);

    mpz_clear(vv);
    return S;
}

Complex zeta_block_mpfr(mpz_t v, unsigned int K, mpfr_t t) {
    //
    // Compute the sum_{k=v}^{v + K - 1) n^{-.5 + it}
    //
    // We use machine doubles for the terms in the sum, but we use
    // mpfr to accurately calculate the quantity t log n mod 2 pi
    // for each term in the sum.
 
    //
    // First we figure out how much precision we will need from mpfr.
    //
    // We want to accurately calculate t log n mod 2pi to 53 bits, which
    // means that we need to compute t log n to 53 + log_2(t log n) bits.
    //
    // For safety we add an extra 2 bits.
    //

    mpfr_t x, y, n, twopi;
    mpfr_init2(x, 53);
    mpfr_init2(y, 53);

    mpfr_log2(x, t, GMP_RNDN);                          // Now x = log2(t)
    mpfr_set_z(y, v, GMP_RNDN);
    mpfr_add_ui(y, y, K, GMP_RNDN);
    mpfr_log(y, y, GMP_RNDN);
    mpfr_log2(y, y, GMP_RNDN);                          // And y = log2(log n) (We calculate these quantities
                                                        // to low precision because we are only interested
                                                        // in their size, really. There is probably 
                                                        // a clever way to do faster using less precision.

    mpfr_add(x, x, y, GMP_RNDN);
    int mod_precision = mpfr_get_ui(x, GMP_RNDN) + 55;  // This is the precision that we need when
                                                        // need to calculate the quantity t log n
                                                        // when we mod by 2 pi

    mpfr_set_z(y, v, GMP_RNDN);
    mpfr_add_ui(y, y, K, GMP_RNDN);
    mpfr_log2(y, y, GMP_RNDN);                          // y = log2(v + K) now.
    int n_precision = mpfr_get_ui(y, GMP_RNDN) + 2;     // This is the precision that we need to exactly
                                                        // represent the largest integer that will occur
                                                        // in the summation.
 
    mpfr_init2(n, n_precision);

    mpfr_init2(twopi, mod_precision);
    mpfr_const_pi(twopi, GMP_RNDN);
    mpfr_mul_ui(twopi, twopi, 2, GMP_RNDN);

    mpfr_clear(x);
    mpfr_init2(x, mod_precision);

    Complex S = 0.0;
    
    mpfr_set_z(n, v, GMP_RNDN);             // The summation starts with n = v
    for(unsigned int k = 0; k <= K-1; k++) {
        mpfr_log(x, n, GMP_RNDN);           // x = log(n)
        mpfr_mul(x, x, t, GMP_RNDN);        // a = t log n
        mpfr_fmod(x, x, twopi, GMP_RNDN);
        Complex z = exp(I * mpfr_get_d(x, GMP_RNDN));
        z = z/sqrt(mpfr_get_d(n, GMP_RNDN));
        S = S + z;
    
        mpfr_add_ui(n, n, 1, GMP_RNDN);     // now on the next iteration, n will be v + k + 1
    }

    mpfr_clear(x);
    mpfr_clear(y);
    mpfr_clear(n);
    mpfr_clear(twopi);

    return S;
}

Complex zeta_block_d_stupid(mpz_t v, int K, mpfr_t t) {
    Double vv = mpz_get_d(v);
    Double tt = mpfr_get_d(t, GMP_RNDN);

    Complex S = 0;
    for(int l = 0; l <= K-1; l++) {
        Double n = vv + l;
        S = S + pow(n, -.5 + I * tt);
    }
    return S;
}

Complex initial_zeta_sum_mpfr(mpz_t M, mpfr_t t) {
    mpfr_t x, y;
//    mpfr_init2(x, 53);
    mpfr_init2(y, 53);
//    mpfr_log2(x, t, GMP_RNDN);
    mpfr_set_z(y, M, GMP_RNDN);
    mpfr_log(y, y, GMP_RNDN);
//    mpfr_log2(y, y, GMP_RNDN);
//    mpfr_add(x, x, y, GMP_RNDN);

    int mod_precision = mpfr_get_exp(t) + mpfr_get_exp(y) + 55;
//    int mod_precision = mpfr_get_ui(x, GMP_RNDN) + 55;  // This is the precision that we need when
                                                        // need to calculate the quantity t log n
                                                        // when we mod by 2 pi

//    mpfr_set_z(y, M, GMP_RNDN);
//    mpfr_log2(y, y, GMP_RNDN);
    int n_precision = mpz_sizeinbase(M, 2) + 2;         // This is the precision that we need to exactly
                                                        // represent the largest integer that will occur
                                                        // in the summation.
    
    
    mpfr_t nn, twopi;
    mpfr_init2(nn, n_precision);
    mpfr_init2(twopi, mod_precision);

    mpfr_const_pi(twopi, GMP_RNDN);
    mpfr_mul_ui(twopi, twopi, 2, GMP_RNDN);

 //   mpfr_clear(x);
    mpfr_init2(x, mod_precision);

    mpz_t n;
    mpz_init(n);

    Complex S1 = 0;

    for(mpz_set_si(n, 1); mpz_cmp(n, M) <= 0; mpz_add_ui(n, n, 1)) {
        mpfr_set_z(nn, n, GMP_RNDN);
        mpfr_log(x, nn, GMP_RNDN);           // x = log(n)
        mpfr_mul(x, x, t, GMP_RNDN);         //  x = t log n

        mpfr_fmod(y, x, twopi, GMP_RNDN);

        Double z = mpfr_get_d(y, GMP_RNDN);
        S1 = S1 + exp(I * z)/sqrt(mpz_get_d(n));

//        unsigned int nn = mpfr_get_ui(n, GMP_RNDN);
//        if(nn % 10000 == 0) {
//            cout << nn << endl;
//        }

//        mpfr_sin_cos(b, a, a, GMP_RNDN);    // a + ib = exp(i t log n)
//        mpfr_sqrt(c, n, GMP_RNDN);          // c = sqrt(n)
//        mpfr_div(a, a, c, GMP_RNDN);
//        mpfr_div(b, b, c, GMP_RNDN);
//        mpfr_add(real_part, real_part, a, GMP_RNDN);
//        mpfr_add(imaginary_part, imaginary_part, b, GMP_RNDN);
    }
    
    mpz_clear(n);
    mpfr_clear(x);
    mpfr_clear(y);
    mpfr_clear(twopi);

    return S1;
}

Complex zeta_sum_mpfr(mpfr_t t) {
    mpfr_t x;
    mpz_t z;

    mpfr_init2(x, mpfr_get_prec(t));
    mpz_init(z);

    mpfr_const_pi(x, GMP_RNDN);                 // x = pi
    mpfr_mul_si(x, x, 2, GMP_RNDN);             // x = 2 pi
    mpfr_div(x, t, x, GMP_RNDN);                // x = t/2pi
    mpfr_sqrt(x, x, GMP_RNDN);                  // x = sqrt(t/2pi)
    mpfr_floor(x, x);                           // x = floor(sqrt(t/2pi))

    mpfr_get_z(z, x, GMP_RNDN);

    Complex S = initial_zeta_sum_mpfr(z, t);

    mpfr_clear(x);
    mpz_clear(z);
    return S;
}


Double C(int n, Double *powers_z);

Double rs_remainder(mpfr_t t) {
    //computes 3 remainder term in riemann-siegel formula
    //can be extended to more than 3 terms easily, if needed
    //everything can be safely done in doubles unless otherwise is stated
    //we're using the notation of the ``amortized complexity..." paper
    //
    //This is adapted from Michael Rubinstein's code for the riemann-siegel formula,
    //the lcalc library
    //

    Double remainderTerms;
    Double c0, c1, c2;
    Double z;
    Double *powers_z;
      
    //variables and calculations below should be done in MPFR
    mpfr_t mp_a, mp_p, twopi;
    mpz_t n1;
    mpfr_init2(mp_a, mpfr_get_prec(t));     // TODO: intelligently choose precision?
    mpfr_init2(twopi, mpfr_get_prec(t));     // TODO: intelligently choose precision?
    mpfr_init2(mp_p, 53);
    mpz_init(n1);

    mpfr_const_pi(twopi, GMP_RNDN);
    mpfr_mul_si(twopi, twopi, 2, GMP_RNDN);

    mpfr_div(mp_a, t, twopi, GMP_RNDN);         // mp_a = t/twopi
    mpfr_sqrt(mp_a, mp_a, GMP_RNDN);            // mp_a = sqrt(t/twopi)

    Double a = mpfr_get_d(mp_a, GMP_RNDN);
    mpfr_get_z(n1, mp_a, GMP_RNDD);                 // n1 = floor(a)

    //a = sqrt(t / (2 * Pi));
    //n1 = floor(a);

    mpfr_frac(mp_p, mp_a, GMP_RNDN);            // mp_p = frac(a)
    //p = a - n1;

    Double p = mpfr_get_d(mp_p, GMP_RNDN);

    //z can be calculated as a double
    z=p-.5;

    //precomputed arrays that probably save time

    powers_z = new Double[51];
    powers_z[0]=1;
    for(int l=1; l<=50; l++)
        powers_z[l] = powers_z[l-1] * z;

    //remainderTerms C(0,p), C(1,p), C(2,p)

    c0=C(0,powers_z);
    c1=C(1,powers_z) * pow(a, -1);
    c2=C(2,powers_z) * pow(a, -2);

    //remainderTerms = C(0,p) + C(1,p)*pow(a,-1) + C(2,p)*pow(a,-2) 

    remainderTerms = c0 + c1 + c2;

    //remainderTerms = remainderTerms * (-1)^(n1+1) / sqrt(a)

    mpz_mod_ui(n1, n1, 2);
    if(mpz_cmp_ui(n1, 1) == 0)
        remainderTerms = remainderTerms * pow(a, -.5);
    else
        remainderTerms = -remainderTerms * pow(a, -.5);

    delete [] powers_z;
    mpz_clear(n1);
    mpfr_clear(twopi);
    mpfr_clear(mp_a);
    mpfr_clear(mp_p);

    return remainderTerms;

}

Double C(int n, Double *powers_z){
    if (n==0) return
        .3826834323650897717284599840304*powers_z[0]
        +1.74896187231008179744118586948533*powers_z[2]
        +2.11802520768549637318456427826417*powers_z[4]
        -.87072166705114807391892407738239*powers_z[6]
        -3.4733112243465167073064116693758*powers_z[8]
        -1.66269473089993244964313630119286*powers_z[10]
        +1.21673128891923213447689352804417*powers_z[12]
        +1.30143041610079757730060538099786*powers_z[14]
        +.03051102182736167242108987123981*powers_z[16]
        -.37558030515450952427981932122934*powers_z[18]
        -.10857844165640659743546975901329*powers_z[20]
        +.05183290299954962337576051067322*powers_z[22]
        +.02999948061990227592040084956912*powers_z[24]
        -.00227593967061256422601994851021*powers_z[26]
        -.00438264741658033830594007013585*powers_z[28]
        -.00040642301837298469930723272116*powers_z[30]
        +.00040060977854221139278910314608*powers_z[32]
        +8.97105799138884129783418195378689e-05*powers_z[34]
        -2.30256500272391071161029452573899e-05*powers_z[36]
        -9.38000660190679248471972940127474e-06*powers_z[38]
        +6.32351494760910750424986123959430e-07*powers_z[40]
        +6.55102281923150166621223123133411e-07*powers_z[42]
        +2.21052374555269725866086890381876e-08*powers_z[44]
        -3.32231617644562883503133517017624e-08*powers_z[46]
        -3.73491098993365608176460476015222e-09*powers_z[48]
        +1.24450670607977391951510000249366e-09*powers_z[50];
    if (n==1) return
        -.05365020525675069405998280791133*powers_z[1]
        +.11027818741081482439896362071917*powers_z[3]
        +1.23172001543152263131956529162206*powers_z[5]
        +1.26349648627994578841755482191213*powers_z[7]
        -1.69510899755950301844944739906596*powers_z[9]
        -2.99987119676501008895548735894141*powers_z[11]
        -.10819944959899208642692257787438*powers_z[13]
        +1.94076629462127126879387632539716*powers_z[15]
        +.78384235615006865328843457488694*powers_z[17]
        -.50548296679003659187902141326173*powers_z[19]
        -.3845072349605797405134273885311*powers_z[21]
        +.03747264646531532067594447494023*powers_z[23]
        +.09092026610973176317258142450576*powers_z[25]
        +.01044923755006450921820113972659*powers_z[27]
        -.01258297965158341649747892224592*powers_z[29]
        -.00339950372115127408505894886137*powers_z[31]
        +.00104109505377148912682954240655*powers_z[33]
        +.00050109490511184868603556526727*powers_z[35]
        -3.95635966900318155954711855696337e-05*powers_z[37]
        -4.76245924535718963865409830268035e-05*powers_z[39]
        -1.85393553380851322734349064569117e-06*powers_z[41]
        +3.19369180800689720404663539343268e-06*powers_z[43]
        +4.09078076085060663265089453677018e-07*powers_z[45]
        -1.54466243325766321284375723273104e-07*powers_z[47]
        -3.46630749176913317222559405934073e-08*powers_z[49];
    if (n==2) return
        .00518854283029316849378458151923*powers_z[0]
        +.00123786335522538984133826974438*powers_z[2]
        -.18137505725166997411491896409414*powers_z[4]
        +.14291492748532126541165603376514*powers_z[6]
        +1.33033917666875653250993329998546*powers_z[8]
        +.35224723534037336775327655505836*powers_z[10]
        -2.4210015958919507237815305433405*powers_z[12]
        -1.67607870225381088533346181492372*powers_z[14]
        +1.36894167233283721842349153807076*powers_z[16]
        +1.55390194302229832214563952655935*powers_z[18]
        -.17221642734729980519582586998918*powers_z[20]
        -.63590680550454309889704902355845*powers_z[22]
        -.09911649873041208105423564341370*powers_z[24]
        +.14033480067387008950738254898316*powers_z[26]
        +.04782352019827292236438803506512*powers_z[28]
        -.01735604064147978079795864709223*powers_z[30]
        -.01022501253402859184447660413126*powers_z[32]
        +.00092741491597948878994270014371*powers_z[34]
        +.00135721943723733853452533619958*powers_z[36]
        +6.41369012029388008996238736394533e-05*powers_z[38]
        -.00012300805698196629883342322937*powers_z[40]
        -1.83135074047892025547675543979621e-05*powers_z[42]
        +7.82162860432262730850139938461872e-06*powers_z[44]
        +2.00875424847599455034985293919157e-06*powers_z[46]
        -3.35327653931857137372749727241453e-07*powers_z[48]
        -1.46160209174182309264510097122760e-07*powers_z[50];
    if (n==3) return
        -.00267943218143891380853967145989*powers_z[1]
        +.02995372109103514963731329491570*powers_z[3]
        -.04257017254182869798501935111688*powers_z[5]
        -.28997965779803887506893209478669*powers_z[7]
        +.48888319992354459725374746407169*powers_z[9]
        +1.23085587639574608119312504336294*powers_z[11]
        -.82975607085274087041796910432976*powers_z[13]
        -2.24976353666656686652045012659903*powers_z[15]
        +.07845139961005471379365473620184*powers_z[17]
        +1.74674928008688940039198666645219*powers_z[19]
        +.45968080979749935109237306173169*powers_z[21]
        -.66193534710397749464339040008983*powers_z[23]
        -.31590441036173634578979632973316*powers_z[25]
        +.12844792545207495988511847476209*powers_z[27]
        +.10073382716626152300969450207513*powers_z[29]
        -.00953018384882526775950465984230*powers_z[31]
        -.01926442168751408889840098069714*powers_z[33]
        -.00124646371587692917124790716458*powers_z[35]
        +.00242439696411030857397215245841*powers_z[37]
        +.00043764769774185701827561290396*powers_z[39]
        -.00020714032687001791275913078304*powers_z[41]
        -6.27434450418651556052610958029804e-05*powers_z[43]
        +1.15753438145956693483789208989316e-05*powers_z[45]
        +5.88385492454037978388597885697078e-06*powers_z[47]
        -3.12467740069633622086961449076033e-07*powers_z[49];
    if (n==4) return
        .00046483389361763381853630462560*powers_z[0]
        -.00402264294613618830391153989145*powers_z[2]
        +.00384717705179612688359130685272*powers_z[4]
        +.06581175135809486002088309200741*powers_z[6]
        -.19604124343694449117695528448205*powers_z[8]
        -.20854053686358853244400012794494*powers_z[10]
        +.95077541851417509458477574151058*powers_z[12]
        +.53415353129148739760517592459894*powers_z[14]
        -1.67634944117634007959116448203404*powers_z[16]
        -1.07674715787512899278784663510432*powers_z[18]
        +1.23533930165659698528788361189251*powers_z[20]
        +1.02578253400572757718348949577914*powers_z[22]
        -.40124095793988544378728137523313*powers_z[24]
        -.50366639951083034479585257591604*powers_z[26]
        +.03573487795502744985807080163387*powers_z[28]
        +.14431763086785416624285239495844*powers_z[30]
        +.01509152741790346941712677290432*powers_z[32]
        -.02609887477919436131761773965448*powers_z[34]
        -.00612662837951926174904909908948*powers_z[36]
        +.00307750312987084118476787782167*powers_z[38]
        +.00115624789340887523161201204220*powers_z[40]
        -.00022775966758472127472807733953*powers_z[42]
        -.00014189637118181444432681579894*powers_z[44]
        +7.46486030795591945312240984450313e-06*powers_z[46]
        +1.24797016454091166174449988846871e-05*powers_z[48]
        +4.86394518400209461907998084746180e-07*powers_z[50];
    if (n==5) return
        .00022686811845737363176557957245*powers_z[1]
        +.00110812468537183880897586725284*powers_z[3]
        -.01621857925555009106408484258686*powers_z[5]
        +.05276503405398741662724126665649*powers_z[7]
        +.02570880200903323999290010111095*powers_z[9]
        -.38058660440806397264435991848146*powers_z[11]
        +.22531987892642315322976926989838*powers_z[13]
        +1.03445733164952217211304499657389*powers_z[15]
        -.55282576970508137898888475296735*powers_z[17]
        -1.52877126410780729962736571714169*powers_z[19]
        +.32828366427719583672031669394059*powers_z[21]
        +1.22911021854008706238425001239677*powers_z[23]
        +.04093693938311529830689289790902*powers_z[25]
        -.55860404726420193442735876775644*powers_z[27]
        -.11241976368059115396788439789609*powers_z[29]
        +.15212677711795591829295940144809*powers_z[31]
        +.05173718845528038784023625510664*powers_z[33]
        -.02561227689700728294043343196050*powers_z[35]
        -.01296367251404617794428713962277*powers_z[37]
        +.00254555748186116327806192744188*powers_z[39]
        +.00211933195108777752885073213414*powers_z[41]
        -9.19139194515677754051761292342159e-05*powers_z[43]
        -.00024413466533855272657049552509*powers_z[45]
        -1.36979826922833871237722416376381e-05*powers_z[47]
        +2.06207850332842380642607669485206e-05*powers_z[49];

    else return 0;
}


Complex rs_rotation(mpfr_t t) {
    //returns rotation factor exp(-i theta(t)) in riemann-siegel formula
    //
    //every thing here is done in doubles unless otherwise is stated
    //
    //

    Double theta;
     
    //first two terms of asymptotic expansion for theta(t)
    //should be done in MPFR, temp should be declated as MPFR variable

    mpfr_t temp, twopi;
    mpfr_init2(twopi, mpfr_get_prec(t));        // TODO: intelligently choose precision?
    mpfr_const_pi(twopi, GMP_RNDN);
    mpfr_mul_si(twopi, twopi, 2, GMP_RNDN);
    
    mpfr_init2(temp, mpfr_get_prec(t));         // TODO: intelligently choose precision?
    
    mpfr_div(temp, t, twopi, GMP_RNDN);         // temp = t/2pi
    mpfr_log(temp, temp, GMP_RNDN);             // temp = log(t/2pi)
    mpfr_sub_ui(temp, temp, 1, GMP_RNDN);       // temp = log(t/2pi) - 1
    mpfr_mul(temp, temp, t, GMP_RNDN);          // temp = t log(t/2pi) - t
    mpfr_div_ui(temp, temp, 2, GMP_RNDN);       // temp = t/2 log(t/2pi) - t/2
    mpfr_fmod(temp, temp, twopi, GMP_RNDN);     // temp = (above) mod 2pi

//    temp =  (t / 2.) * log(t / (2. * pi)) - t / 2.
//    temp = fmod(theta, 2 * pi);

    //a few more terms in asymptotic expansion of theta(t)
    //can be safely done in doubles
    
    Double tt = mpfr_get_d(t, GMP_RNDN);

    theta = mpfr_get_d(temp, GMP_RNDN) - PI/8.0 +  1.0/(tt * 48.0) + 7.0/(tt*tt*tt * 5760.0); 

    //calculate exp(-I*theta(t))

    Complex answer = exp(-I * theta);

    mpfr_clear(twopi);
    mpfr_clear(temp);

    return answer;
}

void compute_Z_from_rs_sum(mpfr_t t0, double delta, int N, Complex * S, Complex * Z) {
    Complex rotation_factor;
    Double remainder_terms;

    mpfr_t t;
    mpfr_init2(t, mpfr_get_prec(t0));
    mpfr_set(t, t0, GMP_RNDN);

    for(int l = 0; l < N; l++) {
        rotation_factor = rs_rotation(t);
        remainder_terms = rs_remainder(t);
        Z[l] = 2 * real(rotation_factor * S[l]) + remainder_terms;
        mpfr_add_d(t, t, delta, GMP_RNDN);
    }
}

void compute_zeta_from_Z(mpfr_t t0, double delta, int N, Complex * Z, Complex * zeta) {
    Complex rotation_factor;
    Double remainder_terms;

    mpfr_t t;
    mpfr_init2(t, mpfr_get_prec(t0));
    mpfr_set(t, t0, GMP_RNDN);

    for(int l = 0; l < N; l++) {
        rotation_factor = rs_rotation(t);
        remainder_terms = rs_remainder(t);
        zeta[l] = 2 * real(rotation_factor * Z[l]) + remainder_terms;
        mpfr_add_d(t, t, delta, GMP_RNDN);
    }

}
