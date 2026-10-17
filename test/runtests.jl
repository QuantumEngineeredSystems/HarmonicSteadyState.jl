using HarmonicSteadyState
using ParallelTestRunner: ParallelTestRunner

# Start with autodiscovered tests
testsuite = ParallelTestRunner.find_tests(@__DIR__)

# Parse arguments
args = ParallelTestRunner.parse_args(ARGS)

if ParallelTestRunner.filter_tests!(testsuite, args)
    delete!(testsuite, "quality/JET")
    # Not part of the suite: MTKHomotopy needs packages outside the test
    # environment, and the Krylov and Float32 files are unfinished.
    delete!(testsuite, "MTKHomotopy")
    delete!(testsuite, "steady_states/krylov")
    delete!(testsuite, "steady_states/float32")
end

ParallelTestRunner.runtests(HarmonicSteadyState, args; testsuite)
