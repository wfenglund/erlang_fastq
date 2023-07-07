# erlang_fastq

In order to run the script, make sure that the fastq-files to be analyzed are
placed in the location '../Raw_data' in relation to the program. In order for
the erlang program to find the fastq-files the forward and reverse reads has
to have names ending in '1.fq.gz' and '2.fq.gz' respectively.
Then navigate to the script and compile and run the program like this:
```bash
$ erlc fastqstats.erl
$ erl -noshell -s fastqstats start -s init stop
```
