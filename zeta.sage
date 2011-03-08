# program to use band-limited function interpolation to compute the zeta
# function using the main sum data precomputed at a grid of points.

from mpmath import siegeltheta, grampoint, mp
import pymongo
from pymongo import Connection

from bisect import bisect
from copy import copy

connection = Connection()
db = connection.zeta
large_values_collection = db.large_values


mp.prec = 300

import sys

R = RealField(300)
C = ComplexField(300)
PI = R(pi)

S = PolynomialRing(R, 'z')
z = S.gen()
C0 = .3826834323650897717284599840304  \
        + 1.74896187231008179744118586948533 * z^2 \
        + 2.11802520768549637318456427826417 * z^4 \
        - .87072166705114807391892407738239 *  z^6 \
        - 3.4733112243465167073064116693758 *  z^8 \
        - 1.66269473089993244964313630119286 * z^10 \
        + 1.21673128891923213447689352804417 * z^12 \
        + 1.30143041610079757730060538099786 * z^14 \
        + .03051102182736167242108987123981 *  z^16 \
        - .37558030515450952427981932122934 *  z^18 \
        - .10857844165640659743546975901329 *  z^20 \
        + .05183290299954962337576051067322 *  z^22 \
        + .02999948061990227592040084956912 *  z^24 \
        - .00227593967061256422601994851021 *  z^26 \
        - .00438264741658033830594007013585 *  z^28 \
        - .00040642301837298469930723272116 *  z^30 \
        + .00040060977854221139278910314608 *  z^32 \
        + 8.97105799138884129783418195378689e-05 * z^34 \
        - 2.30256500272391071161029452573899e-05 * z^36 \
        - 9.38000660190679248471972940127474e-06 * z^38 \
        + 6.32351494760910750424986123959430e-07 * z^40 \
        + 6.55102281923150166621223123133411e-07 * z^42 \
        + 2.21052374555269725866086890381876e-08 * z^44 \
        - 3.32231617644562883503133517017624e-08 * z^46 \
        - 3.73491098993365608176460476015222e-09 * z^48 \
        + 1.24450670607977391951510000249366e-09 * z^50

C1 =    - .05365020525675069405998280791133  * z        \
        + .11027818741081482439896362071917  * z^3      \
        + 1.23172001543152263131956529162206 * z^5      \
        + 1.26349648627994578841755482191213 * z^7      \
        - 1.69510899755950301844944739906596 * z^9      \
        - 2.99987119676501008895548735894141 * z^11     \
        - .10819944959899208642692257787438  * z^13      \
        + 1.94076629462127126879387632539716 * z^15      \
        + .78384235615006865328843457488694  * z^17      \
        - .50548296679003659187902141326173  * z^19      \
        - .3845072349605797405134273885311   * z^21      \
        + .03747264646531532067594447494023  * z^23      \
        + .09092026610973176317258142450576  * z^25      \
        + .01044923755006450921820113972659  * z^27      \
        - .01258297965158341649747892224592  * z^29      \
        - .00339950372115127408505894886137  * z^31      \
        + .00104109505377148912682954240655  * z^33      \
        + .00050109490511184868603556526727  * z^35      \
        - 3.95635966900318155954711855696337e-05  * z^37      \
        - 4.76245924535718963865409830268035e-05  * z^39      \
        - 1.85393553380851322734349064569117e-06  * z^41      \
        + 3.19369180800689720404663539343268e-06  * z^43      \
        + 4.09078076085060663265089453677018e-07  * z^45      \
        - 1.54466243325766321284375723273104e-07  * z^47      \
        - 3.46630749176913317222559405934073e-08  * z^49

C2 =      .00518854283029316849378458151923       * z^0       \
        + .00123786335522538984133826974438       * z^2       \
        - .18137505725166997411491896409414       * z^4       \
        + .14291492748532126541165603376514       * z^6       \
        + 1.33033917666875653250993329998546      * z^8       \
        + .35224723534037336775327655505836       * z^10      \
        - 2.4210015958919507237815305433405       * z^12      \
        - 1.67607870225381088533346181492372      * z^14      \
        + 1.36894167233283721842349153807076      * z^16      \
        + 1.55390194302229832214563952655935      * z^18      \
        - .17221642734729980519582586998918       * z^20      \
        - .63590680550454309889704902355845       * z^22      \
        - .09911649873041208105423564341370       * z^24      \
        + .14033480067387008950738254898316       * z^26      \
        + .04782352019827292236438803506512       * z^28      \
        - .01735604064147978079795864709223       * z^30      \
        - .01022501253402859184447660413126       * z^32      \
        + .00092741491597948878994270014371       * z^34      \
        + .00135721943723733853452533619958       * z^36      \
        + 6.41369012029388008996238736394533e-05  * z^38      \
        - .00012300805698196629883342322937       * z^40      \
        - 1.83135074047892025547675543979621e-05  * z^42      \
        + 7.82162860432262730850139938461872e-06  * z^44      \
        +2.00875424847599455034985293919157e-06   * z^46      \
        -3.35327653931857137372749727241453e-07   * z^48      \
        -1.46160209174182309264510097122760e-07   * z^50


def rotation_factor(t):
    theta = t/2 * log(t/(2 * PI)) - t/2 - PI/8 + 1/(48 * t) + 7/(t^3 * 5760)
    return exp(-I * theta)

def remainder_terms(t):
    a = sqrt(t/(2 * PI))
    N = floor(a)
    p = a - N - 1/2

    return_value = C0(p) + C1(p)/a + C2(p)/a^2

    if is_odd(N):
        return_value = return_value/sqrt(a)
    else:
        return_value = -return_value/sqrt(a)

    return return_value

def N_approx(t):
    return RealField(200)(siegeltheta(t)/pi + 1)
    pass

class ZetaData:
    def __init__(self, data, type = "file"):
        if type == "file":
            input_file = open(data, 'r')
            self.t0 = R(input_file.readline().strip())
            self.delta = R(input_file.readline().strip())
            self.data = []
            next = input_file.readline().strip()
            while(len(next) > 0):
                a, b = next[1:-1].split(',')
                self.data.append(C(a,b))
                next = input_file.readline().strip()
            input_file.close()
            self.tmax = self.t0 + self.delta * (len(self.data) - 1)
        elif type == "large_value_entry":
            db_entry = large_values_collection.find_one( {'t' : data} )
            self.t0 = R(db_entry['t'])
            self.delta = R(db_entry['delta'])
            self.data = []
            for x,y in db_entry['values']:
                self.data.append( C(x,y) )

            self.tmax = self.t0 + self.delta * (len(self.data) - 1)
        elif type == "dict":
            self.t0 = data['t0']
            self.delta = data['delta']
            self.data = copy(data['data'])
            self.tmax = self.t0 + self.delta * (len(self.data) - 1)

    def update_db_entry(self, collection):
        t_string = self.t0.str(no_sci = 2, skip_zeroes = True)
        if t_string[-1] == '.':
            t_string = t_string[:-1]
        current_entry = collection.find_one( {'t' : t_string} )
        max_value = float(max( (abs(y) for (x,y) in self.Z_values_from_array() )))
        values = [ (float(real(x)), float(imag(x))) for x in self.data ]
        if current_entry:
            ID = current_entry['_id']
            collection.save( { '_id' : ID, 't' : t_string, 'delta' : float(self.delta), 'max_value' : max_value , 'values' : values} )
        else:
            collection.save( { 't' : t_string, 'delta' : float(self.delta), 'max_value' : max_value, 'values' : values} )


    def Z_values_from_array(self):
        Z_value_pairs = []
        for n in range(len(self.data)):
            t = self.t0 + n * self.delta
            Z_value = 2 * real(self.data[n] * rotation_factor(t)) + remainder_terms(t)
            Z_value_pairs.append( (n * self.delta, Z_value) )
        return Z_value_pairs

    def __call__(self, t):
        return self.Z_value(t)

    def Z_value(self, t, offset=True):
        RR = ComplexField(300)
        CC = ComplexField(300)
        I = CC.0
        PI = RR(pi)
        t = RR(t)
        t0 = RR(self.t0)

        if not offset:
            t = t - t0

        delta = RR(self.delta)
        data_range = self.delta * (len(self.data) - 1)
        if t < 0 or t > data_range:
            raise ValueError

        #if(n0 == 0):
        #    return self.data[0] + t * (self.data[1] - self.data[0])/self.delta
        #if(n0 == len(self.data) - 1):
        #    N = len(self.data) - 2
        #    d = t - N * self.delta
        #    return self.data[N] + d * (self.data[N + 1] - self.data[N])/self.delta

        if data_range - t > t:
            interpolation_range = srange(0, 2 * floor(t/delta))
        else:
            interpolation_range = srange(2 * ceil(t/delta) - len(self.data), len(self.data))

        N = len(interpolation_range)
        if N > 40:
            interpolation_range = interpolation_range[N/2 - 20:N/2 + 20]

        tau = 1/2 * log(self.t0/(2 * PI))
        beta = PI/delta
        Lambda = (beta + tau)/2
        epsilon_1 = (beta - tau)/2

        c = RR(len(interpolation_range))/2 * PI * epsilon_1/beta

        alpha = tau

        S = 0
        for n in interpolation_range:
            u = n * delta - t
            z = RR(self.data[n]) * exp(-I * alpha * u) * sinc(Lambda * u) * blfi_kernel(u, c, epsilon_1)
            S = S + z

        S = S * Lambda/beta

        X = rotation_factor(t0 + t)

        return RealField(53)(2 * real(S * X) + remainder_terms(t0 + t))

    def zeta_value(self, t):
        RR = RealField(200)
        t0 = RR(self.t0)
        t = RR(t)
        S = RR(self.Z_value(t))
        return ComplexField(53)( S * rotation_factor(t0 + t) )

    def zeta_value_separated(self, t):
        # same as above, but returns the real and imaginary values separated in a tuple
        RR = RealField(200)
        t0 = RR(self.t0)
        t = RR(t)
        S = RR(self.Z_value(t))
        zeta_value = ComplexField(53)( S * rotation_factor(t0 + t) )
        return real(zeta_value), imag(zeta_value)

    def make_animation(self, start, end, delta, location, bound):
        frame_number = 1

        t = start
        L = []
        while(t < end):
            L.append(self.zeta_value_separated(t))
            P = list_plot(L, plotjoined=True)
            filename = location + ("/frame%010d" % frame_number) + ".png"
            P.save(filename, figsize=[10,10], xmin=-bound, xmax=bound, ymin=-bound, ymax=bound)
            frame_number += 1
            t += delta

    def calculate_N(self, t):
        # we choose a gram point g near t and then try to use turing's method
        # to prove that S(g) = 0, so that N(g) = N_approx(g)
        

        RR = RealField(200)
        t = self.t0 + RR(t)

        N = floor(N_approx(t))
        g = RR(grampoint(N))

        delta = RR(.001)

        # we need (-1)^N Z(g) > 0, so we keep advancing to the next gram point if this is not the case
        
        while (-1)^N * self.Z_value(g, offset=False) <= .001:
            N = N + 1
            g = grampoint(N)

        starting_N = N

        print "using N = ", N, "; g = ", g, " as starting gram point."

        g_list = [RR(g)]
        h_list = [RR(0)]

        upper_bound_points = [RR(g)]

        a = RR(2.067) #
        b = RR(.059)  # values from Trudgian's paper

        S_upper_bound = 100

        # we know that S(g) is an even integer. we start by trying to bound it from above by something less than 2

        while S_upper_bound > 1.95:
            N = N + 1
            g2 = RR(grampoint(N))
            if g_list[-1] + h_list[-1] < g2 and (-1)^N * self.Z_value(g2, offset=False) > 0:
                g_list.append(g2)
                h_list.append(0)
                S_upper_bound = RR(1 + (a + b * log(g2) + sum(h_list))/(g2 - g))
                print "Found upper bound of ", S_upper_bound
                upper_bound_points.append(g2)
            else:
                h2 = g_list[-1] + h_list[-1] - g2 + delta
                while (-1)^N * self.Z_value(g2 + h2, offset = False) > 0:
                    h2 = h2 + delta
                    print "trying h =", h2
                
                print "appending point", N, ",", h2

                g_list.append(g2)
                h_list.append(h2)

                upper_bound_points.append(g2 + h2)

                print "Not computing upper bound because h is not zero"

        print "Successfully proved that S(g) <= 0"

        # now we do the same thing moving to the left, and trying to prove a lower bound for S(g)

        S_lower_bound = -100
        N = starting_N
        lower_bound_points = [g]
        g_list = [g]
        h_list = [0]
        while S_lower_bound < -1.95:
            N = N - 1
            g2 = RR(grampoint(N))
            if g_list[-1] + h_list[-1] > g2 and (-1)^N * self.Z_value(g2, offset=False) > 0:
                g_list.append(g2)
                h_list.append(0)
                S_lower_bound = RR(-1 - (a + b * log(g2) - sum(h_list))/(g - g2))
                print "Found lower bound of ", S_lower_bound
                lower_bound_points.append(g2)
            else:
                h2 = g_list[-1] + h_list[-1] - g2 - delta
                while (-1)^N * self.Z_value(g2 + h2, offset = False) > 0:
                    h2 = h2 - delta
                    print "trying h =", h2
                
                print "appending point", N, ",", h2

                g_list.append(g2)
                h_list.append(h2)

                print "Not computing lower bound because h is not zero"

                lower_bound_points.append(g2 + h2)

        print "Successfully proved that S(g) >= 0. Huzzah! "

        return (starting_N, RR(g - self.t0), upper_bound_points, lower_bound_points)

    def find_zeros(self, start, end, delta):
        t1 = start
        t2 = t1 + delta
        zeros = []
        count = 0
        while(t2 <= end):
            if self.Z_value(t1) * self.Z_value(t2) < 0:
                z = find_root(self, t1, t2)
                zeros.append(z)
                count = count + 1
                percent_done = RR((t2 - start)/(end - start) * 100)
                print "Found zero number", count, "at", z, "; ", percent_done, "percent done"
            t1 = t2
            t2 = t1 + delta

        return zeros

    def create_S_t(self, start, starting_index, zeros):
        def S(t):
            if t < start:
                return 0
            num_zeros = starting_index + bisect(zeros, t) + 1
            return num_zeros - RealField(200)(N_approx(self.t0 + RealField(200)(t)))
            
        #return plot(S, start, end, rgbcolor=(1, 0, 0), plot_points=5000)
        return S


def plot_max_values():
    entries = large_values_collection.find()
    L = []
    log10 = RR(10).log()
    for entry in entries:
        t = RR(entry['t'])
        m = RR(entry['max_value'])
        L.append( point((t.log(), m.log())))

    return sum(L)



def update_all_entries(collection):
    t_list = []
    entries = collection.find(fields = ['t'])
    for entry in entries:
        t_list.append(entry['t'])
    
    for t in t_list:
        print "loading data for t =", t
        Z = ZetaData(t, type = "large_value_entry")
        print "updating entry for t =", t
        Z.update_db_entry(collection)

def update_all_entries_missing_max(collection):
    t_list = []
    entries = collection.find(fields = ['t', 'max_value'])
    for entry in entries:
        if not 'max_value' in entry:
            t_list.append(entry['t'])

    for t in t_list:
        print "loading data for t =", t
        Z = ZetaData(t, type = "large_value_entry")
        print "updating entry for t =", t
        Z.update_db_entry(collection)

def blfi_kernel(u, c, epsilon_1):
    x = c^2 - (epsilon_1 * u)^2
    if abs(x) < .00001:
        w = 1/39916800*x^5 + 1/362880*x^4 + 1/5040*x^3 + 1/120*x^2 + 1/6*x + 1
    else:
        w = sinh(sqrt(x))/sqrt(x)
    z = w * c / cosh(c)
    return z

def sinc(x):
    if abs(x) < .00001:
        return 1/120*x^4 - 1/6*x^2 + 1
    else:
        return sin(x)/x

def make_pictures():
    data_directory = 'data2/'
    output_directory = 'pictures2/'
    #exps = [18, 19, 20, 21, 22, 23, 24, 25, 26, 27]
    exps = [28,]

    for n in exps:
        datafile = data_directory + 'zeta_1e' + str(n) + '.data'
        filename_base = output_directory + 'z_1e' + str(n)
        large_picture_name = filename_base + '_large.png'
        small_picture_name = filename_base + '_small.png'
        Z = ZetaData(datafile)
        L = [(x-2, y) for (x,y) in Z.Z_values_from_array()]
        
        P = list_plot(L, plotjoined = True)

        print 'Writing image to', large_picture_name
        sys.stdout.flush()
        P.save(large_picture_name, figsize=[60,20])

        print 'Writing image to', small_picture_name
        sys.stdout.flush()
        P.save(small_picture_name, figsize=[10,5])

def make_pictures_from_database(collection):
    output_location = "static/images"
    
    t_list = []
    entries = collection.find(fields = ['t'])
    for entry in entries:
        t_list.append(entry['t'])
    
    for t in t_list:
        print "Processing data for t =", t
        large_location = os.path.join( output_location, "large", t + ".png")
        small_location = os.path.join( output_location, "small", t + ".png")
        axis_zoom_location = os.path.join( output_location, "axis_zoom", t + ".png")

        Z = ZetaData(t, type="large_value_entry")

        plot_begin = .1
        plot_end = Z.tmax - Z.t0 - .1

        plot_created = False
        if not os.path.exists(small_location):
            if not plot_created:
                P = plot(lambda t : Z(t), plot_begin, plot_end)
                plot_created = True
            P.save(small_location, figsize=[20,10])

        if not os.path.exists(large_location):
            if not plot_created:
                P = plot(lambda t : Z(t), plot_begin, plot_end)
                plot_created = True
            P.save(large_location, figsize=[100,30])

        if not os.path.exists(axis_zoom_location):
            if not plot_created:
                P = plot(lambda t : Z(t), plot_begin, plot_end)
                plot_created = True
            P.save(axis_zoom_location, figsize=[100,30], ymax=10, ymin=-10)

def process_incoming():
    data_directory = 'data7/'
    output_directory = 'pictures3/'
    filenames = os.listdir(data_directory)

    for input_file in filenames:
        datafile = data_directory + input_file
        Z = ZetaData(datafile)
        Z.update_db_entry(large_values_collection)
        t = ZZ(floor(Z.t0))
        #L = [(x-2, y) for (x,y) in Z.Z_values_from_array()]
        #L = Z.Z_values_from_array()
        #P = list_plot(L, plotjoined = True)
        print 'calculating plot for t = ', t
        sys.stdout.flush()
        P = plot(lambda t : Z(t), .1, 39.9)
        
        large_picture_name = output_directory + 'z_' + str(t) + '_large.png'
        vertical_zoom_picture_name = output_directory + 'z_' + str(t) + '_vertical_zoom.png'
        axis_zoom_picture_name = output_directory  + 'z_' + str(t) + '_axis_zoom.png'
        small_picture_name = output_directory + 'z_' + str(t) + '_small.png'
        medium_picture_name = output_directory + 'z_' + str(t) + '_medium.png'

        print 'Writing image to', large_picture_name
        sys.stdout.flush()
        P.save(large_picture_name, figsize=[30,10], dpi=200)

        print 'Writing image to', small_picture_name
        sys.stdout.flush()
        P.save(small_picture_name, figsize=[5,2.5], dpi=200)

        #print 'Writing image to', vertical_zoom_picture_name
        #sys.stdout.flush()
        #P.save(vertical_zoom_picture_name, figsize=[20,300])

        print 'Writing image to', axis_zoom_picture_name
        sys.stdout.flush()
        P.save(axis_zoom_picture_name, figsize=[30,10], ymin=-5, ymax=5, dpi=200)

        print 'Writing image to', medium_picture_name
        sys.stdout.flush()
        P.save(medium_picture_name, figsize=[8,4], ymin=-5, ymax=5, dpi=200)

def get_max_values():
    data_directory = 'data/all/'
    filenames = os.listdir(data_directory)

    max_values = []

    for input_file in filenames:
        datafile = data_directory + input_file
        print "Processing file", input_file
        sys.stdout.flush()
        Z = ZetaData(datafile)
        t = ZZ(Z.t0)

        L = [abs(y) for (x, y) in Z.Z_values_from_array()]
        max_values.append((t, max(L)))

    return max_values



def find_candidate_large_value():
    #primes = [2,3,5,7,11,13,17]
    P = prime_range(600)
    N = len(P)
    m = 84 
    r = 15 
    possible_t = []

    for n in range(70, 90):
        primes = P[:n]
        weights = [p^(1/4) for p in primes]
        for m in [70, 71, 72, 73, 74, 75, 76, 77, 78]:
            for r in [11, 12, 13, 14, 15, 16, 17, 18, 19]:
                A = copy(MatrixSpace(ZZ, n+1).zero())
                for k in xrange(0, n):
                    A[0,k] = floor(weights[k] * RR(primes[k]).log() * 2^(m - r))

                A[0,n] = 1

                for k in xrange(0, n):
                    A[k + 1, k] = floor(2 * pi * weights[k] * 2^m)

                #A = permutation_action(SymmetricGroup(n).random_element(), A)

                B = A.LLL(delta=.9)
                #print B.str()
                R = RealField(200)
                for k in xrange(0, n + 1):
                    possible_t.append(R(B[n, k]/2^r))

    return possible_t

def find_candidates_for_large_value(repeat=1):
    possible_t = []
    X = var('X')
    
    p_start = 110 
    p_end = 150 

    euler_product1 = 1
    euler_product2 = 1
    for p in prime_range(nth_prime(p_start) + 1):
        print p
        euler_product1 = euler_product1 / (1 - 1/p^(1/2 + i * X))

    euler_product2 = euler_product1

    for p in prime_range(nth_prime(p_start) + 1, nth_prime(p_end) + 1):
        euler_product2 = euler_product2 / (1 - 1/p^(1/2 + i * X))

    euler_product1 = fast_callable(euler_product1, domain=ComplexField(150), vars='X')
    euler_product2 = fast_callable(euler_product2, domain=ComplexField(150), vars='X')

    for l in xrange(repeat):
        n = ZZ.random_element(p_start, p_end)
        m = ZZ.random_element(90, 130)
        r = ZZ.random_element(40, 50)
        delta = RR.random_element(.9, .9999)

        last_prime = nth_prime(n)
        primes = prime_range(last_prime + 1)

        weights = [p^(-1/4) for p in primes]
        #print weights

        A = copy(MatrixSpace(ZZ, n+1).zero())
        for k in xrange(0, n):
            A[0,k] = floor(weights[k] * primes[k].log() * 2^(m - r))

        A[0,n] = 1

        for k in xrange(0, n):
            A[k + 1, k] = floor(2 * pi * weights[k] * 2^m)

        B = A.LLL(delta=delta)

        R = RealField(200)
        for k in xrange(0, n + 1):
            t = R(abs(B[n,k])/2^r)
            v1 = euler_product1(t)
            v2 = euler_product2(t)
            possible_t.append( (abs(v2), abs(v1), t, floor(abs(t)) - 2) )
        
        possible_t.sort()
        possible_t = possible_t[-500:]
        print "On try number", l, " best candidates so far are"
        for _ in possible_t[-10:]:
            print _
        sys.stdout.flush()

    return possible_t


def find_candidate_large_value2():
    #primes = [2,3,5,7,11,13,17]
    primes = prime_range(30)

    weights = [p^(1/4) for p in primes]

    m = 82 
    r = 15 
    n = len(primes)

    A = copy(MatrixSpace(ZZ, n+1).zero())
    for k in xrange(0, n):
        A[k,0] = floor(weights[k] * primes[k] * 2^(m - r))

    A[n,0] = 1

    for k in xrange(0, n):
        A[k, k+1] = floor(2 * pi * weights[k] * 2^m)

    B = A.LLL()


    possible_t = []
    R = RealField(200)
    for k in xrange(0, n + 1):
        possible_t.append(R(B[k, n]/2^r))

    return possible_t

def euler_product(X, domain=CC):
    p = 2
    S = 1
    x = var('x')
    while p < X:
        S = S/(1 - 1/p^(1/2 + i * x))
        p = next_prime(p)

    f = fast_callable(S, vars=[x], domain=domain)
    return f


def euler_product_values(t_list, X):
    p = 2
    S = 1

    x = var('x')
    while p < X:
        S = S/(1 - 1/p^(1/2 + i * x))
        #S = S * CC((1 - p^(-i * t)/sqrt(p))^(-1))
        p = next_prime(p)

    f = fast_callable(S, vars=[x], domain=CC)
    L = [ (abs(f(t)), t) for t in t_list ]
    return L

def blah():
    L = find_candidate_large_value()
    M = euler_product_values(L, 500)

    #N = [y for (x,y) in M]
    M.sort()
    
    for x, y in M:
        print y, x

def blah2(x):
    return floor(x) - 2

def blah3(location):
    work_units = os.listdir(location)
    work_units.sort()

    print work_units

    results_created = False

    filecount = 0

    for filename in work_units:
        print "processing", os.path.join(location, filename)
        w = open(os.path.join(location, filename), 'r')
        first_line = w.readline()
        t_str, start_str, length_str, N_str, delta_str, filename_str, cputime_str, realtime_str = first_line.split()

        t0 = RealField(300)(t_str)
        N = Integer(N_str)
        delta = R(delta_str)
        if not results_created:
            results = [C(0)] * N
            results_created = True

        count = 0
        for line in w:
            x, y = line.strip()[1:-1].split(',')
            x = R(x)
            y = R(y)
            #results[count] += CC(line)
            results[count] += C( (x,y) )
            count = count + 1

        w.close()

        filecount = filecount + 1
        if filecount > 253:
        #if True:
            Z = ZetaData({"t0" : t0, "delta" : delta, "data" : results}, type="dict")
            L = Z.Z_values_from_array()
            P = list_plot(L, plotjoined = True)
            #P.save("partial_sums/%08d.png" % filecount, ymax=80, ymin=-80, figsize=[30,20])
            P.save("partial_sums/%08d.png" % filecount, figsize=[10,5])

    Z = ZetaData({"t0" : t0, "delta" : delta, "data" : results}, type="dict")
    #L = Z.Z_values_from_array()
    #P = list_plot(L, plotjoined = True)
    P = plot(lambda t : Z(t), 1, 39)
    P.save("partial_sums/final.png", ymax=15, ymin=-15, figsize=[30,20])


def examine_euler_product(t, X):
    C = ComplexField(200)
    I = C("(0, 1)")
    euler_product = C(1)
    n = 0
    prev = 1
    for p in prime_range(nth_prime(X) + 1):
        n = n + 1
        prev = euler_product
        euler_product = euler_product / (1 - 1/p^(1/2 + I * t))
        average = (prev * euler_product)^(1/2)
        if p % 5000 == 3:
            print n, p, CC(euler_product), average

    print p, CC(euler_product), average

def remove_duplicates(zero_list, delta):
    output = [zero_list[0]]
    for z in zero_list[1:]:
        if abs(z - output[-1]) > delta:
            output.append(z)

    return output

def compute_spacings(zeros):
    s = [ (y - x) for (x,y) in zip(zeros[:-1], zeros[1:]) ]
    s.sort()
    return s

def partition_zeros(zeros, g1, g2):
    z1 = [z for z in zeros if z < g1]
    z2 = [z for z in zeros if z > g1 and z < g2]
    z3 = [z for z in zeros if z > g2]

    return z1, z2, z3

def make_plot_for_msri():
    t0 = RealField(300)("9178358656494989336431259004785")
    offset = 20.337429013215523
    values = load("zeta_values")
    
    abs_values = [(x - offset, abs(y)) for x, y in values]
    P = list_plot(abs_values, plotjoined=True)
    P.show(figsize=[30,15], xmin=-2, xmax=2)
