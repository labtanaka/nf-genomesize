#!/usr/bin/env nextflow

/*
* MIT License
*
* Copyright (c) 2019 Sergej Nowoshilow (sergej.nowoshilow@imp.ac.at)
*
* Permission is hereby granted, free of charge, to any person obtaining a copy
* of this software and associated documentation files (the "Software"), to deal
* in the Software without restriction, including without limitation the rights
* to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
* copies of the Software, and to permit persons to whom the Software is
* furnished to do so, subject to the following conditions:
*
* The above copyright notice and this permission notice shall be included in all
* copies or substantial portions of the Software.
*
* THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
* IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
* FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
* AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
* LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
* OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
* SOFTWARE.
*/

params.reads = '/groups/tanaka/Data/sequences/Illumina/Genome_correction/2nd_round/84393'
params.chunkSize = 10_000_000
params.canonical = 'yes'
params.kmer = 25

def canonical = (params.canonical && params.canonical.toLowerCase() == 'yes') ? '--canonical' : ''


Channel
  .fromPath(params.reads + '/**.fastq*')
  .splitFastq(by: params.chunkSize, compress: false, file: true)
  .set { ch_fastq_files }


process countKMers {

  input:
    file(fastqfile) from ch_fastq_files

  output:
    file('counts.jf') into ch_counts

  script:
  """
  jellyfish count --size=15G --mer-len=${params.kmer} --threads=${task.cpus} --output=counts.jf ${canonical} --out-counter-len=2 ${fastqfile}
  """
}


ch_counts
  .collect()
  .collectFile { files -> ["countslist", files.collect{ it.toString() }.join('\n')] }
  .set { ch_counts_list }

process mergeCounts {

  publishDir "./results"

  input:
    file(list) from ch_counts_list

  output:
    file('merged.counts.jf') into ch_merged_counts

  script:
  """
  jellyfish merge --output=merged.counts.jf \$(cat ${list} | tr "\n" " ")
  """ 
}


process buildHistogram {

  publishDir "./results"

  input:
    file(counts) from ch_merged_counts

  output:
    file('counts.hist') into ch_counts_hist

  script:
  """
  jellyfish histo --threads=${task.cpus} --output=counts.hist ${counts}
  """
}