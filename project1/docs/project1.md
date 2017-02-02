# Xiaoya Zhou's miRNA mircroarray data analysis

<!-- TOC START min:1 max:3 link:true update:true -->
- [Xiaoya Zhou's miRNA mircroarray data analysis](#xiaoya-zhous-mirna-mircroarray-data-analysis)
  - [Count numbers for each type of miRNAs/probes](#count-numbers-for-each-type-of-mirnasprobes)
  - [Extract expression matrix from raw data](#extract-expression-matrix-from-raw-data)
  - [Normalize expression matrix](#normalize-expression-matrix)

<!-- TOC END -->

## Count numbers for each type of miRNAs/probes

```perl
perl -F, -lane'!($.>= 12 and $F[3] =~ /miR|let/) and next; $count{(split/[-_]/, $F[3])[1]}++ }{ print qq{$_ | $count{$_}} for sort {$count{$b} <=> $count{$a}} keys %count' 'Raw Intensity File.csv'
```

STDOUT:

Type | Counts
---- | ------
miR | 8200
miRPlus | 72
let | 72

## Extract expression matrix from raw data

> Perl-oneliner I used below removed headers (comment lines at the start of the file) and abundant columns.

> Known that each kind of probe has four replication on the microarray, I chose the probe set with **the highest variance** for each miRNA as representation. Replicated miRNAs were merged by **median** of each sample by *KangChen Bio-tech*, which seems like a self-evidenct method, but has fatal limitation. See [**Question: Handling Duplicate Probe Expression Values In Spotted Cdna Microarray** on **Biostars**](https://www.biostars.org/p/51756/#51875) for details.

```perl
perl -MStatistics::Descriptive -F, -lane'!($.>= 12 and $F[3] =~ /miR|let/) and next; @tmp = split/,/, $_; $stat->add_data(@tmp[23..$#tmp]); $tvar = $stat->variance(); if ($tvar > $var{$F[3]}) {$var{$F[3]} = $tvar; $line{$F[3]} = join(qq{\t}, @tmp[3,23..$#tmp])} $stat->clear()}{ BEGIN{$stat = Statistics::Descriptive::Sparse->new()} print $line{$_} for keys %line' 'Raw Intensity File.csv' > expr
```

Through checking the results, file `expr` has 2081 lines while the `Raw Intensity File.csv` file has 2085 probe records. What's the **four** ones remaining? Comparison applied to show difference.

```bash
# To get a sorted miRNA name list
cut -f1 expr | sort > 1
cut -d, -f4 'Raw Intensity File.csv' | sort > 2
# Then I remove headers manually (actually it should be deleted by `sed`:) )
diff 1 2 > 3
cat 3
```

The result of `diff`:

```pre
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
```

These miRNAs had more than one one probes. Considering the requirement of identifying differentially expressed miRNAs and this rare scenario, my procedure may works well up to now.

## Normalize expression matrix

To illustrate raw data distribution, a density plot and a boxplot were generated.

```R
setwd('~/github.com/bioinformatist/research_projects/project1/')
library(data.table)
library(cowplot)

DT.expr1 <- fread('expr')
sample.names <- c('miRNAs', '37', '40', '58', '59', '63', '69', '49', '50', '64', '67', '68', '79', '70', '72', '73', '75', '90', '91')
setnames(DT.expr1, sample.names)
DT.expr1 <- melt(DT.expr1, id = 1, variable.name = 'sample.name')
```

> Here, both package `reshape` and `reshape2` have function `melt`, and `data.table` has its own `melt.data.table` function alias `melt`. How can we deal with them?

```R
# To specify the package that you want to use, the syntax is: packagename::functionname()
data.table::melt
reshape2::melt
# If you always want to use function in one, you can define your own function as follows:
melt <- data.table::melt
```

> And keep in mind that the order of loading the packages makes a difference, i.e. the package that gets loaded last will mask the functions in packages loaded earlier.

Let's move on. Make graphics now:

```R
density.raw <- ggdraw(ggplot(data = DT.expr1, aes(value, colour = sample.name)) + geom_density() + scale_x_continuous(trans = 'log2')) + draw_label("Draft for \n Peng's Lab!", angle = 45, size = 80, alpha = .2)
save_plot('figures/density.raw.png', density.raw, base_height = 8.5, base_width = 11)
boxplot.raw <- ggdraw(ggplot(data = DT.expr1, aes(sample.name, value)) + geom_boxplot(notch = TRUE) + scale_y_continuous(trans = 'log2')) + draw_label("Draft for \n Peng's Lab!", angle = 45, size = 80, alpha = .2)
save_plot('figures/boxplot.raw.png', density.raw, base_height = 8.5, base_width = 11)
```




Why use such size for figures?

> About figure size: Each figure should be able to fit on a single 8.5 x 11 inch page. Please do not send figure panels as individual files. We use three standard widths for figures: 1 column, 85 mm; 1.5 column, 114 mm; and 2 column, 174 mm (the full width of the page). Although your figure size may be reduced in the print journal, please keep these widths in mind. For Previews and other three-column formats, these widths are also applicable, though the width of a single column will be 55 mm. --From [Cell Press Digital Image Guidelines (click to see details)](http://www.cell.com/figureguidelines).
