


"""
examples

easy
    report_data/multiqc.log
    report_data/multiqc_.+\.txt
    (?P&lt;name&gt;.+)\.png
    mqc_(?P&lt;designation&gt;.+)\.txt

medium
    split_bam(?P&lt;designation&gt;.+)\.bam
    (?P&lt;designation&gt;.+)\.report\.tsv
    (?P&lt;designation&gt;.+)\.sdf$
    (?P&lt;designation&gt;.+)_ss\.ps
    (?P&lt;designation&gt;.+_vs_.+)
    gffcmp\.(?P&lt;designation&gt;.+)\.tmap
    .*?\.(?P&lt;designation&gt;.*)\.fasta
    .*?\.(?P&lt;designation&gt;.*)\.cons\.tax\.summary
    .*?\.column\.dist\.(?P&lt;designation&gt;.*)\.temp

hard
    (?P&lt;identifier_0&gt;\S+\ \S+)\ single 1\.(?P&lt;ext&gt;.*)
    (?P&lt;identifier_0&gt;.*?) index (0|1)\.(?P&lt;ext&gt;.*)
    (?P&lt;identifier_0&gt;\S+ \S+) (single (0|1)|(forward|reverse) 0)\.(?P&lt;ext&gt;.*)

"""