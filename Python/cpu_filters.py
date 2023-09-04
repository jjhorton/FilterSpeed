import os
import sys
import time
import json
import math
import array
import tempfile
import random
import collections

from gnuradio import gr
from gnuradio import audio, analog, digital, filter, blocks

#import clenabled
import time

filter_length =100

################################################################################

# Benchmark parameters

# Duration of each benchmark trial
BENCH_TRIAL_DURATION    = 5.0
# Number of benchmark trials to average
BENCH_NUM_TRIALS        = 11
# Benchmark Suite
BenchmarkSuite = []

################################################################################

# Decorator for defining benchmarks in the suite
def benchmark(test_name, block_name):
    def wrapped(f):
        BenchmarkSuite.append((test_name, block_name, f))
        return f
    return wrapped

################################################################################

@benchmark("FIR Filter (Complex taps, Complex input", "filter.fir_filter_ccc")
def test_fir_filter_ccc():
    top = gr.top_block()
    src = blocks.null_source(gr.sizeof_gr_complex)
    firfilter = filter.fir_filter_ccc(1, [complex(random.random(), random.random()) for _ in range(filter_length)])
    probe = blocks.probe_rate(gr.sizeof_gr_complex)
    top.connect(src, firfilter, probe)

    return top, probe

@benchmark("FIR Filter (1 Thread, FFT, Complex taps, Complex input", "filter.fft_filter_ccc")
def test_fft_filter_ccc():
    top = gr.top_block()
    src = blocks.null_source(gr.sizeof_gr_complex)
    firfilter = filter.fft_filter_ccc(1, [complex(random.random(), random.random()) for _ in range(filter_length)],nthreads=1)
    probe = blocks.probe_rate(gr.sizeof_gr_complex)
    top.connect(src, firfilter, probe)

    return top, probe

@benchmark("FIR Filter (2 Thread, FFT, Complex taps, Complex input", "filter.fft_filter_ccc")
def test_fft_filter_ccc():
    top = gr.top_block()
    src = blocks.null_source(gr.sizeof_gr_complex)
    firfilter = filter.fft_filter_ccc(1, [complex(random.random(), random.random()) for _ in range(filter_length)],nthreads=2)
    probe = blocks.probe_rate(gr.sizeof_gr_complex)
    top.connect(src, firfilter, probe)

    return top, probe

@benchmark("FIR Filter (4 Thread, FFT, Complex taps, Complex input", "filter.fft_filter_ccc")
def test_fft_filter_ccc():
    top = gr.top_block()
    src = blocks.null_source(gr.sizeof_gr_complex)
    firfilter = filter.fft_filter_ccc(1, [complex(random.random(), random.random()) for _ in range(filter_length)],nthreads=4)
    probe = blocks.probe_rate(gr.sizeof_gr_complex)
    top.connect(src, firfilter, probe)

    return top, probe

################################################################################

# Benchmark runner

if __name__ == '__main__':
    # If a test name was specified, filter the benchmark suite by fuzzy-matching
    # by test name
    if len(sys.argv) > 1:
        MatchedBenchmarkSuite = []

        for benchmark in BenchmarkSuite:
            if benchmark[0].lower().find(sys.argv[1].lower()) >= 0:
                MatchedBenchmarkSuite.append(benchmark)

        BenchmarkSuite = MatchedBenchmarkSuite

    benchmark_results = {
        'version': gr.version(),
        'parameters': {
            'num_trials': BENCH_NUM_TRIALS,
            'trial_duration': BENCH_TRIAL_DURATION
        },
        'benchmarks': []
    }

    for index, benchmark in enumerate(BenchmarkSuite):
        for filter_length in range( 1, 1002, 50):
            test_name, block_name, test_factory = benchmark
            test_name = test_name + (", length %i)" % filter_length)
            sys.stderr.write("Running benchmark {}/{} \"{}\"\n".format(index+1, len(BenchmarkSuite), test_name))

            samples_per_second, bytes_per_second = [], []

            # Run each trial
            for trial in range(BENCH_NUM_TRIALS):
                # Create the test top block
                test_top, test_probe = test_factory()

                # Run the trial
                test_top.start()
                time.sleep(BENCH_TRIAL_DURATION)
                test_top.stop()

                trial_samples_per_second = test_probe.rate()
                trial_bytes_per_second = trial_samples_per_second * test_probe.input_signature().sizeof_stream_item(0)

                sys.stderr.write("\tTrial {} - {:.3f} MS/s, {:.2f} MiB/s\n".format(trial+1, trial_samples_per_second/1e6, trial_bytes_per_second/1048576))

                samples_per_second.append(trial_samples_per_second)
                bytes_per_second.append(trial_bytes_per_second)

                time.sleep(1)

            # Compute means
            mean_samples_per_second = sum(samples_per_second)/BENCH_NUM_TRIALS
            mean_bytes_per_second = sum(bytes_per_second)/BENCH_NUM_TRIALS

            # Compute standard deviations
            stdev_samples_per_second = math.sqrt(sum([(e - mean_samples_per_second)**2 for e in samples_per_second])/BENCH_NUM_TRIALS)
            stdev_bytes_per_second = math.sqrt(sum([(e - mean_bytes_per_second)**2 for e in bytes_per_second])/BENCH_NUM_TRIALS)

            sys.stderr.write("\tAverage - {:.3f} MS/s, {:.1f} MiB/s\n".format(mean_samples_per_second/1e6, mean_bytes_per_second/1048576))
            sys.stderr.write("\t  Stdev - {:.3f} MS/s, {:.1f} MiB/s\n".format(stdev_samples_per_second/1e6, stdev_bytes_per_second/1048576))

            # Add it to our table
            benchmark_results['benchmarks'].append({
                'name': test_name,
                'block_name': block_name,
                'results': {
                    'samples_per_second': mean_samples_per_second,
                    'samples_per_second_stdev': stdev_samples_per_second,
                    'bytes_per_second': mean_bytes_per_second
                }
            })

    print(json.dumps(benchmark_results))
