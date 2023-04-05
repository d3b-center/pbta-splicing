
function stdev_mean(TITLE) { return(DATA[TITLE,"T"]/DATA[TITLE,"C"]); }
function stdev_count(TITLE) { return(DATA[TITLE,"C"]-1); }

# Run first for data 1-n to get means
function stdev_pass1(TITLE, VAL) {
        DATA[TITLE,"T"] += VAL+0;
        DATA[TITLE,"C"] ++;
}

# Run second for data 1-n to get standard deviations
function stdev_pass2(TITLE,VAL,X) {
        X = stdev_mean(TITLE) - VAL+0;
        DATA[TITLE,"D"] += X*X;
}

# Final result after both passes
function stdev(TITLE) { return(sqrt(DATA[TITLE,"D"] / stdev_count(TITLE)));}

