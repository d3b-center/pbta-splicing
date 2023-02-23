# Integrate titles and skip first line
FNR==1 {        for(N=2; N<=NF; N++) COL[N]=$N ; MAX=NF; next   }

# First pass, calculate means and skip to next line
NR==FNR {       for(N=2; N<=NF; N++) stdev_pass1(N, $N); next }

# Second pass, means are now valid, calculate deviation and correlation
{
        for(N=2; N<=NF; N++) stdev_pass2(N, $N);
        for(N=2; N<=NF; N++) for(M=N; M<=NF; M++)
                CORR[N,M]+=(stdev_mean(N) - ($N+0)) * (stdev_mean(M) - ($M+0));
}

END { # Print final data
        for(N=2; N<=MAX; N++)   for(M=N; M<=MAX; M++)
        print COL[N], COL[M], CORR[N,M] / (stdev_count(N)*stdev(N)*stdev(M));
}

# Not a typo, awk is fed inputfile twice.
# This avoids needing to store the entire massive file in memory.

