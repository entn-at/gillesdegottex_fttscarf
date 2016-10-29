#include <iostream>
#include <iterator>
#include <deque>
#include <fstream>
#include <sstream>
using namespace std;

#include <fftscarf.h>
using namespace fftscarf;

#include <boost/chrono/system_clocks.hpp>
#include <boost/random/normal_distribution.hpp>

#include <boost/random/mersenne_twister.hpp>
static boost::random::mt19937 rndgen(std::time(0));

#include "stream.h"

template<typename FFTPlanType>
static void benchmark(const string& name, int minlN=2, int maxlN=17){
    try{

    std::cout << "Testing " << name << " ..." << std::endl;
    std::cout << "Impl. type N MFlops dur5[s] dur50[s] dur95[s] acc5[s] acc50[s] acc95[s] nbruns" << endl;

//    fftscarf::write_compile_info(cout);
    boost::mt19937 rnd_engine;
    boost::random::normal_distribution<FFFLOAT> rnd_normal_distrib;
    boost::chrono::system_clock::time_point tstart;
    boost::chrono::system_clock::time_point tend;
    std::vector<std::complex<FFFLOAT> > spec;
    std::vector<std::complex<FFFLOAT> > spec_ref;
    std::vector<FFFLOAT> inframe, outframe;

//     fftscarf::FFTPlanImplementationFFTW3 fft_ref(true); // TODO
    fftscarf::FFTPlanImplementationOoura fft_ref(true); // TODO
    FFTPlanType fft(true);
    FFTPlanType ifft(false);
    bool spec_verify = false; // TODO

    std::deque<double> durations; // [mus]
    std::deque<double> accuracies;

    ofstream resultfile("benchmark.log");

    // Need to start at 32 for ffts
    vector<int> Ns;
    for(int lN=minlN; lN<=maxlN; ++lN)
        Ns.push_back(std::pow(2, lN));
    for(int Ni=0; Ni<Ns.size(); ++Ni){
        int N = Ns[Ni];

        durations.clear();
        accuracies.clear();
        
        for(size_t runi=0; runi<10000; ++runi) {
            // Window the frame
            inframe.resize(N);   // The input frame
            for(size_t n=0; n<N; ++n)
                inframe[n] = rnd_normal_distrib(rnd_engine);

            if(spec_verify)
                fft_ref.dft(inframe, spec_ref, N);

            tstart = boost::chrono::system_clock::now();

            fft.dft(inframe, spec, N);
            ifft.idft(spec, outframe, N);

            tend = boost::chrono::system_clock::now();

            boost::chrono::nanoseconds nanosec = tend-tstart;
            durations.push_back(1e-9*nanosec.count()/2); // per DFT

//             outframe.resize(inframe.size()); // TODO Something weired happens with FFTReal

            if(spec_verify){
                // Verify: sig->spec == specref
                FFFLOAT sn = 0.0;
                for(size_t i=0; i<inframe.size(); ++i)
                    sn += abs(spec_ref[i]-spec[i])*abs(spec_ref[i]-spec[i]);
                sn = sqrt(sn/inframe.size());

    //                if(sn>1000*fftscarf::eps){
    //                    DCOUT << sn << " eps=" << 1000*fftscarf::eps << std::endl;
    //                    DCOUT << "m_spec_ref=" << spec_ref << endl;
    //                    DCOUT << "m_spec=" << spec << endl;
    //                }
    //                assert(sn<1000*fftscarf::eps);
            }

            // Verify: sig->spec->sig' == sig
            // Using relative RMS
            FFFLOAT sqerr = 0.0;
            FFFLOAT sqin = 0.0;
            for(size_t i=0; i<inframe.size(); ++i){
                sqerr += (inframe[i]-outframe[i])*(inframe[i]-outframe[i]);
                sqin += inframe[i]*inframe[i];
            }
            
            accuracies.push_back(sqrt(sqerr/sqin));
        }

        std::sort(durations.begin(), durations.end());
        double dur_per5 = durations[int(0.05*durations.size())];
        double dur_per50 = durations[int(0.50*durations.size())];
        double dur_per95 = durations[int(0.95*durations.size())];

        std::sort(accuracies.begin(), accuracies.end());
        double acc_per5 = accuracies[int(0.05*accuracies.size())];
        double acc_per50 = accuracies[int(0.50*accuracies.size())];
        double acc_per95 = accuracies[int(0.95*accuracies.size())];

        stringstream result;
        result << name << " sr " << N << " " << 0.5*5*N*log2(N)/(1e6*dur_per50) << " " << dur_per5 << " " << dur_per50 << " " << dur_per95 << " " << acc_per5 << " " << acc_per50 << " " << acc_per95 << " " << durations.size() << endl;

        // Output result in both file and standard output
        resultfile << result.str() << flush;
        cout << result.str() << flush;
    }

    }catch(string err){
        std::cout << err << std::endl;
    }
}