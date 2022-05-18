
# Calculate % of achieved memory bandwidth
To display the below inline math correctly, make sure to have the [xhub](https://github.com/nschloe/xhub) chrome extension as github does not support LaTex Math yet.

The formula for calculating Bandwidth of the 256x256 input is as follows:

```math
Bandwidth = \frac{size\_of\_one\_grid * (timestep(params.maxIters))}{Elapsed\ Compute\ Time}
```

```math
    = \frac{256*256 * 9 (speeds) * 4 (bytes\_per\_float) * 2(cells+tmp\_cells) * 2(read+write) * 80000(timestep)}{35.7s}
```

```math
    = \frac{2.36(MiBytes) * 4 * 80000}{35.7s} = 21.25GB/s
```
 The maximum achievable L2 memory bandwidth is 84.88 GB/s, from here we can calculate the fraction of STREAM bandwidth = $`\frac{21.25}{84.88} = 25\%`$.

# Calculate parallel efficiency
The formula for calculating parallel efficiency is
```math
Parallel\ Efficiency (PE) = \frac{s}{n}
```
(where $`s`$ is ration of speedup, $`n`$ is the processor/core count).

Take the result from [MPI_speedup.svg](../plot/MPI/svg/MPI_speedup.svg) as an exaxmple, calculating PE of a 1024x1024 grid running with MPI and 112 cores:
```math
PE_{mpi} = \frac{40}{112} = 35.71\%
```
