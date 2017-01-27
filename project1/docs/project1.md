# Xiaoya Zhou's miRNA mircroarray data analysis

## Count numbers for each type of miRNAs/probes

```perl
perl -F, -lane'!($.>= 12 and $F[3] =~ /miR|let/) and next; $count{(split/[-_]/, $F[3])[1]}++ }{ print qq{$_ | $count{$_}} for sort {$count{$b} <=> $count{$a}} keys %count' 'Raw Intensity File.csv'
```

STDOUT:

miR | 8200
miRPlus | 72
let | 72

## Extract expression matrix from raw data

> Perl-oneliner I used below removed headers (comment lines at the start of the file) and abundant columns.
> Known that each kind of probe has four replication on the microarray, I chose the probe set with the highest variance for each miRNA as representation.

```perl
perl -MStatistics::Descriptive -F, -lane'!($.>= 12 and $F[3] =~ /miR|let/) and next; @tmp = split/,/, $_; $stat->add_data(@tmp[23..$#tmp]); $tvar = $stat->variance(); if ($tvar > $var{$F[3]}) {$var{$F[3]} = $tvar; $line{$F[3]} = join(qq{\t}, @tmp[3,23..$#tmp])} $stat->clear()}{ BEGIN{$stat = Statistics::Descriptive::Sparse->new()} print $line{$_} for keys %line' 'Raw Intensity File.csv' > expr
```

Through checking the results, file `expr` has 2081 lines while the `Raw Intensity File.csv` file has 2085 probe records. What's the **four** ones remaining? I did some simple statistics.

```bash
cut -f1 expr | sort > 1
cut -d, -f4 'Raw Intensity File.csv' | sort > 2
# Then I remove headers manually (actually it should be deleted by `sed`:) )
diff 1 2 > 3
cat 3
```

The result of `diff`:

0a1
>
256a258
> hsa-miR-134-3p
685a688
> hsa-miR-328-5p
860a864
> hsa-miR-381-5p
1931a1936
> hsa-miR-874-5p

These miRNAs had more than one one probes. Considering the requirement of identifying differentially expressed miRNAs and this rare scenario, my procedure may works well up to now.
