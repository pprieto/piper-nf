/*
 * Copyright (c) 2013, Centre for Genomic Regulation (CRG) and the authors.
 *
 *   This file is part of 'Piper-NF'.
 *
 *   Piper-NF is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 *
 *   Piper-NF is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with Piper-NF.  If not, see <http://www.gnu.org/licenses/>.
 */


import groovyx.gpars.dataflow.operator.DataflowEventAdapter
import groovyx.gpars.dataflow.operator.DataflowProcessor
import com.google.common.collect.Multiset
import com.google.common.collect.HashMultiset
import nextflow.util.CacheHelper
import java.nio.file.Files


/* 
 * Main Piper-NF pipeline script
 *
 * @authors
 * Giovanni Bussotti <giovannibussotti@gmail.com>
 * Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 * Pablo Prieto <poena.funesta@gmail.com> 
 *
 *
 * Defines the pipeline parameters.
 * The values in the 'params' map can be overridden on the command line by specifying a
 * option prefixed with a double '-' char, for example
 *
 * $ nextflow piper.nf --query=<path to your file name>
 *
 */

params.queryChunkSize = 100
params.query = 'tutorial/5_RNA_queries.fa'
params.genomesDb = 'db'
params.resultDir = 'result'
params.blastStrategy = 'ncbi-blast'     // the blast tool to be used, choose between: ncbi-blast, wu-blast
params.alignStrategy = 'slow_pair'      // defines the T-Coffee alignment method
params.exonerateSuccess = '1'
params.exonerateMode = 'exhaustive'
params.exonerateChunkSize = 200
params.referenceGtf = null

// these parameters are mutually exclusive
// Input genome can be specified by
// - genomes-file: a file containing the list of genomes FASTA to be processed
// - genomes-list: a comma separated list of genomes FASTA file
// - genomes-folder: a directory containing a folder for each genome FASTA file
//params['genomes-file'] = null
params['genomes-list'] = null
params['genomes-folder'] = "tutorial/genomes/"

queryFile   = file(params.query)
dbPath      = file(params.genomesDb)
alnPath     = file(params.resultDir+"/aln")
gtfPath     = file(params.resultDir+"/gtf")
fastaPath   = file(params.resultDir+"/fasta")
mfaPath     = file(params.resultDir+"/mfa")

if( !dbPath.exists() ) {
    log.warn "Creating genomes-db path: $dbPath"
    if( !dbPath.mkdirs() ) {
        exit 1, "Cannot create genomes-db path: $dbPath -- check file system permissions"
    }
}

log.info "P I P E R - RNA mapping pipeline - ver 1.1"
log.info "=========================================="
log.info "query               : ${queryFile}"
log.info "genomes-db          : ${dbPath}"
log.info "query-chunk-size    : ${params.queryChunkSize}"
log.info "result-dir          : ${params.resultDir}"
log.info "blast-strategy      : ${params.blastStrategy}"
log.info "align-strategy      : ${params.alignStrategy}"
log.info "exonerate-success:  : ${params.exonerateSuccess}"
log.info "exonerate-mode:     : ${params.exonerateMode}"
log.info "exonerate-chunk-size: ${params.exonerateChunkSize}"
log.info "pool-size           : ${config.poolSize}"
log.info "\n"

/* Define path to reference Gtf */
if( params['referenceGtf']) {
    refAbsPath = file(params.referenceGtf)
}
/*
 * Find out all the genomes files in the specified directory.
 *
 * More in detail teh 'sourceGenomesPath' points to a directory having a
 * sub-folder for each genome it is required to process.
 *
 * Each sub-folder must contain the genome FASTA file to be processed.
 *
 * The sub-folder name is used to identify the genome in the computation.
 *
 * All the genomes names found in this path are put in a list named 'formatName',
 * which control the pipeline execution.
 *
 */

allGenomes = [:]

// when the provided source path is a FILE
// each line represent the path to a genome file
if( params['genomes-file'] ) {
    def genomesFile = file(params['genomes-file'])
    if( genomesFile.empty() ) {
        exit 1, "Not a valid input genomes descriptor file: ${genomesFile}"
    }

    allGenomes = parseGenomesFile(dbPath, genomesFile, params.blastStrategy, params.referenceGtf)
}

else if( params['genomes-list'] ) {
   allGenomes = parseGenomesList(dbPath, params['genomes-list'], params.blastStrategy)
}

else if( params['genomes-folder'] ) {
    def sourcePath = file(params['genomes-folder'])
    if( !sourcePath.exists() || sourcePath.empty() ) {
        exit 4, "Not a valid input genomes folder: ${sourcePath}"
    }

    allGenomes = parseGenomesFolder(dbPath, sourcePath, params.blastStrategy)
}

else {
    exit 5, "No input genome(s) provided -- Use one of the following CLI options 'genomes-file' or 'genomes-list' or 'genomes-folder' "
}

if( !allGenomes ) {
    exit 6, "No genomes found in path"
}

allGenomes.each { name, entry ->
    log.info "Validating genome: $name -- file: ${entry.genome_fa}"
    if( !entry.genome_fa.exists() ) {
        exit 3, "Missing genome file: ${entry.genome_fa}"
    }
}


/*
 * Split the query input file in many small files (chunks).
 *
 * The number of sequences in each chunk is controlled by the parameter 'queryChunkSize'
 * The chunk files are saved in a local folder define by the variable 'querySplits'
 *
 */

// create a folder that may be cached, using the 'queryFile' and the number chunks as cache key
querySplits = cacheableDir([queryFile, params.queryChunkSize])

if( querySplits.empty() ) {
    log.info "Splitting query file: $queryFile .."
    chunkCount=0
    queryFile.chunkFasta( params.queryChunkSize ) { sequence ->
        def file = querySplits.resolve("seq_${chunkCount++}")
        file.text = sequence
    }
    log.info "Created $chunkCount input chunks to path: ${querySplits}"
}
else {
    log.info "Cached query splits > ${querySplits.list().size()} input query chunks"
}


// use a set since there should be not repetition
allQueryIDs = new HashSet()

Path queryEntries = cacheableDir(queryFile)

queryFile.chunkFasta() { String chunk ->
    String queryId = chunk.readLines()[0].replaceAll( /^>(\S*).*$/, '$1' )

    allQueryIDs << queryId
    // store the chunk to a file named as the 'queryId'
    def fileEntry = queryEntries.resolve(queryId)
    if( fileEntry.empty() ) {
        fileEntry.text = chunk
    }
}


/*
 * Create the required databases (BLAST,CHR) if they does not exists.
 *
 * This task is executed for each genome in the list 'formatName'
 * The tasks 'sends' out the name of the genome to be processed
 * by the next step in the pipeline using the variable 'blastName'
 */


def sed_cmd = (System.properties['os.name'] == 'Mac OS X' ? 'gsed' : 'sed')
def split_cmd = (System.properties['os.name'] == 'Mac OS X' ? 'gcsplit' : 'csplit')

process format {
    input :
    val formatName using allGenomes.keySet()
    output:
    val formatName using blastName
    
    """
    set -e
    NAME=${formatName}
    FASTA=${allGenomes[formatName].genome_fa}
    CHR_DB=${allGenomes[formatName].chr_db}
    BLAST_DB=${allGenomes[formatName].blast_db}

    ## Create the BLAST db if they does not exist
    if [[ ! `ls -A ${BLAST_DB} 2>/dev/null` ]]; then

        ## Create the target folder
        mkdir -p ${BLAST_DB}

        ## Format the BLAST DB
        x-format.sh ${params.blastStrategy} ${FASTA} ${BLAST_DB}
    fi


    ## Create the CHR database if does not exist
    if [[ ! `ls -A ${CHR_DB} 2>/dev/null` ]]; then

        ## split the fasta in a file for each sequence 'seq_*'
        ${split_cmd} ${FASTA} '%^>%' '/^>/' '{*}' -f seq_ -n 5

        ## create the target folder
        mkdir -p ${CHR_DB}

        ## rename and move to the target folder
        for x in seq_*; do
        SEQID=`grep -E "^>" $x | ${sed_cmd} -r 's/^>(\\S*).*/\\1/'`
        mv $x ${CHR_DB}/$SEQID;
        done

    fi
    echo $NAME > blastName
    """
}


/*
 * Implements the BLAST step
 */

process blast {
    input:
    each blastId using blastName
    file blastQuery using (querySplits.listFiles())

    output:
    val blastId using exonerateId
    file blastQuery using exonerateQuery
    file '*.mf2' using blastResult
    
    """
    x-blast.sh '${params.blastStrategy}' ${allGenomes[blastId].blast_db} ${blastQuery} > ${blastId}.mf2
    """

}

/*
 * == Blast post-process
 *
 * Split blastResult to small chunks chunks containing at most 'exonerateChunkSize' lines,
 * this chunks feed the exonerate step
 */

exonerate_in = channel()
operator( inputs: [exonerateId, exonerateQuery, blastResult], outputs: [exonerate_in] ) { specieId, fileQuery, fileBlast ->

    fileBlast.chunkLines( size: params.exonerateChunkSize, autoClose: false  ) { lines ->
        // create the chunk file
        def fileChunk = cacheableFile( lines, 'chunk' )
        if( !fileChunk.exists() ) {
            fileChunk.text = lines
        }

        // create 3-tuple to feed to 'exonerate' step
        def id = specieId.trim()
        exonerate_in << [ specie: id, query: fileQuery, chunk: fileChunk, chr_db: allGenomes[id].chr_db ]
    }

}


/*
 * Collect the BLAST output chunks and apply the 'exonerate' function
 */

process exonerate {
    input:
    val exonerate_in

    output:
    file '*.fa' using exonerateOut
    file '*.gtf' using exonerateGtf
    
    """
    specie='${exonerate_in.specie}'
    chr=${exonerate_in.chr_db}
    ## apply exonerate
    exonerateRemapping.pl -query ${exonerate_in.query} -mf2 ${exonerate_in.chunk} -targetGenomeFolder \$chr -exonerate_lines_mode ${params.exonerateMode} -exonerate_success_mode ${params.exonerateMode} -ner no

    mv chunk.fa \${specie}.fa
    mv chunk.ex.gtf \${specie}.ex.gtf

    """
}


/*
 * post-process 'exonerate' result
 */
homologOut      = channel()
normalizedFasta = channel()
normalizedGtf = channel()
normalizationDone = val()
dir = tempDir()
log.debug "Folder exonerateHits: ${dir}"


def foo() {
    normalizationDone << 1
}


def listener = new DataflowEventAdapter() {
    @Override
    public void afterStop(final DataflowProcessor processor) {
        foo()
    }

    public boolean onException(DataflowProcessor processor, java.lang.Throwable e) {
        e.printStackTrace()
        return true
    }
}

Multiset hitSet = HashMultiset.create()

operator( inputs:[exonerateOut, exonerateGtf], outputs: [normalizedFasta, normalizedGtf, homologOut], maxForks: 1, listeners: [listener] ) { fasta, gtf ->

    def specie = fasta.baseName
    def replace = []
    def newFasta = cacheableFile( fasta )

    fasta.chunkFasta(autoClose:false) { seq ->
 
        // parse the sequence id
        def seqId = seq.readLines()[0]
        def matcher = (seqId =~ />(.*)_(hit\d*)(.*)/ )
        def (queryId, hitName, extra) = matcher[0][1..3]
        def sequence = seq.readLines()[1..-1].join('\n')

        // create a multi-fasta file for each 'queryId'
        if( !allQueryIDs.contains(queryId) ) {
            println "Skipping queryId: $queryId -- since it's not contained in the source query"
        }

        log.debug "Processing queryId: ${queryId}"
        def file = dir.resolve("${queryId}.mfa")
        log.debug "MFa: ${file}"
        if( !file.exists() ) {
            // the very fist time prepend the sequence in the query file
            log.debug "File created: ${file}"
            file = Files.createFile(file)
            file.text = queryEntries.resolve(queryId).text
            // note: the file is bound over the channel here, to be sure
            // to send it out just one time
            normalizedFasta << file
        }

        // update the hit name
        def key = [specie, queryId]
        def count = hitSet.add(key, 1) +1
        def newHit = "hit$count"
        if( hitName != newHit ) {
            log.debug "Replacing hitName: $hitName with: $newHit using key: $key"
            replace << [queryId: queryId, oldHit: hitName, newHit: newHit ]
            hitName = newHit
        }

        // now append the query content
        file << ">${queryId}_${hitName}${extra}_${specie}\n"
        file << sequence
        file << '\n'

        newFasta << ">${queryId}_${hitName}${extra}\n"
        newFasta << sequence
        newFasta << '\n'

        homologOut << newFasta
    }

    // normalizing hitNames
    if( replace ) {
        def str = gtf.text
        replace?.each {
            log.debug "Replacing hitName: $it in GTF file: $gtf"
            def pattern = "hitName \"${it.queryId}_${it.oldHit}\";"
            str = str.replaceAll( ~/$pattern/, "hitName \"${it.queryId}_${it.newHit}\";" )
        }

        def newGtf = cacheableFile( gtf )
        newGtf.text = str
        log.debug "Updated GTF file: $newGtf"

        gtf = newGtf
    }
    log.debug "End of operator"
    // send out the 'gtf' file
    normalizedGtf << gtf
}

homolog2Blast   = channel()

homologOut.each { fastaFile ->
    if ( fastaFile.size() == 0 ) return
    def name = fastaFile.name
    def targetFile = fastaPath.resolve(name)
    targetFile << fastaFile.text
    homolog2Blast << fastaFile
}

process orthologs {
    input:
    file homolog2Blast
    output:
    file '*.orthologs.fa' using orthologOut

    """
    set -e
    specie='${homolog2Blast.baseName}'
    touch ${specie}.orthologs.fa

    if [[ ${params['referenceGtf']} != 'null' ]]
    then
        reference_gtf=${refAbsPath}
        ## split the fasta in a file for each sequence 'seq_*'
        ${split_cmd} $homolog2Blast '%^>%' '/^>/' '{*}' -f seq_ -n 5 -s

        ## rename and move to the target folder
        for x in seq_*; do

            x-ortholog.sh '${params.blastStrategy}' $dbPath/reference/${params.blastStrategy}-db $x > bestResult
            cp bestResult best$x

            QUERY=`grep '>' $x | ${sed_cmd} -r 's/^>(.*)_hit\\d*.*\$/\\1/'`
            grep "\\"$QUERY\\"" $reference_gtf > query_gtf
            cp query_gtf gtf$x
            overlap.pl bestResult query_gtf | grep -w 1  | sort | uniq > overlap
            cp overlap overlap$x
            check=`cat overlap`
            if [[ $check -eq 1 ]]
            then
                cat $x >>${specie}.orthologs.fa
            fi
        done
    fi
    """
}

/* Copy the fasta of the 1-to-1 orthologs */
resultDir = file(params.resultDir)
resultDir.with {
    if( !empty() ) { deleteDir() }
    mkdirs()
}

if( !fastaPath.mkdirs() ) {
    exit 1, "Cannot create alignments path: $fastaPath -- check file system permissions"
}

process align {
    input:
    val normalizationDone
    file normalizedFasta

    output:
    file '*.aln' using alignment

    """
    t_coffee -method ${params.alignStrategy} -in $normalizedFasta -n_core 1
    """
}

/*
 * Copy alignment files (MFA and ALN)
 */

if( !alnPath.mkdirs() ) {
    exit 1, "Cannot create alignments path: $alnPath -- check file system permissions"
}
if ( !mfaPath.mkdirs()) {
    exit 1, "Cannot create MFA path: $mfaPath -- check file system premissions"
}

alignme = channel()

alignment.each{ alnFile ->
    def mfaFile = dir.resolve(alnFile.baseName+".mfa")
    if ( mfaFile.size() == 0 ) return

    def targetMfa = mfaPath.resolve(mfaFile.name)
    targetMfa << mfaFile.text

    if ( alnFile.size() == 0 ) return

    def targetAln = alnPath.resolve(alnFile.name)
    targetAln << alnFile.text
    alignme << alnFile
}

process similarity(merge:true) {
    input:
    file alignme

    output:
    file '*' using similarity joint true

    """
    t_coffee -other_pg seq_reformat -in $alignme -output sim > ${alignme.baseName}
    """
}


/* Copy the GFT files produces by the Exonerate steps into the result (current) folder */
if( !gtfPath.mkdirs() ) {
    exit 1, "Cannot create alignments path: $alnPath -- check file system permissions"
}

normalizedGtf.each { sourceFile ->
    if( sourceFile.size() == 0 ) return

    def name = sourceFile.name
    def targetFile = gtfPath.resolve(name)
    targetFile << sourceFile.text
}

orthologOut.each { sourceFile ->
    if ( sourceFile.size() == 0 ) return
    def name = sourceFile.name
    def targetFile = fastaPath.resolve(name)
    targetFile << sourceFile.text
}

/*
 * Compute the similarity Matrix
 */
process matrix {
    echo true
    input:
    file similarity

    output:
    file simMatrix

    """
    echo '\n====== Pipe-R sim matrix ======='
    mkdir data
    mv ${similarity} data
    sim2matrix.pl -query $queryFile -data_dir data -genomes_dir $dbPath | tee simMatrix
    echo '\n'
    """
}

simMatrixFile = simMatrix.val
simMatrixFile.copyTo( resultDir.resolve('simMatrix.csv') )


// ----==== utility methods ====----

def parseGenomesFile(def dbPath, def sourcePath, String blastStrategy, String referenceGtf) {

    def result = [:]

    // parse the genomes input file files (genome-id, path to genome file)
    int count=0
    sourcePath.eachLine { line ->

        def genomeId
        def path

        def items = line.trim().split(/\s+/)
        if( items.size() > 1 ) {
            count++
            (path, genomeId) = items
        }
        else if( items.size() ==1 && items[0] ){
            count++
            genomeId = "gen${count}"
            path = items[0]
        }
        else {
            return
        }
        if( (genomeId == 'reference' && referenceGtf != null) || (genomeId != 'reference')) {
                result[ genomeId ] = [
                        genome_fa: file(path),
                        chr_db: dbPath.resolve("${genomeId}/chr"),
                        blast_db: dbPath.resolve("${genomeId}/${blastStrategy}-db")
                ]
        }
    }

    result
}



def parseGenomesList(def dbPath, String genomesList, String blastStrategy) {

    def count=0
    def files = genomesList.split(',').collect { file(it.trim()) }
    def result = [:]

    files.each { genomeFile ->

        def genomeId = "gen${++count}"
        result[ genomeId ] = [
                genome_fa: genomeFile,
                chr_db: dbPath.resolve("${genomeId}/chr"),
                blast_db: dbPath.resolve("${genomeId}/${blastStrategy}-db")
            ]
    }
    result
}

def parseGenomesFolder(def dbPath, def sourcePath, String blastStrategy) {

    def result = [:]

    sourcePath.eachDir { Path path ->
        def fasta = path.listFiles().find { Path file -> file.name.endsWith('.fa') }
        if( fasta ) {
            result[ path.name ] = [
                    genome_fa: fasta,
                    chr_db: absPath.resolve("${path.name}/chr"),
                    blast_db: absPath.resolve("${path.name}/${blastStrategy}-db")
                ]
        }
    }
    result
}

// ----===== TEST ====-------

def void testParseGenomesFile() {

    def db = file('db')
    def source = file('test-source')
    try {
        source.text =
            '''
            x/file1.fa
            y/file2.fa   genx
            z/file3.fa
            '''

        def result = parseGenomesFile(db, source, 'wu-blast')

        assert result.size() == 3

        assert result['gen1'].genome_fa == file('x/file1.fa')
        assert result['genx'].genome_fa == file('y/file2.fa')
        assert result['gen3'].genome_fa == file('z/file3.fa')

        assert result['gen1'].chr_db == file(db, 'gen1/chr')
        assert result['gen3'].chr_db == file(db, 'gen3/chr')
        assert result['genx'].chr_db == file(db, 'genx/chr')

        assert result['gen1'].blast_db == file(db, 'gen1/wu-blast-db')
        assert result['genx'].blast_db == file(db, 'genx/wu-blast-db')
        assert result['gen3'].blast_db == file(db, 'gen3/wu-blast-db')

    }
    finally {
        source.delete()
    }
}



def void testParseGenomesList() {

      def db = file('db')

      // call the function to test
      def result = parseGenomesList(db, 'alpha.fa, beta.fa, delta.fa', 'x-blast')

      // verify result
      assert result.size() == 3

      assert result['gen1'].genome_fa == file('alpha.fa')
      assert result['gen2'].genome_fa == file('beta.fa')
      assert result['gen3'].genome_fa == file('delta.fa')

      assert result['gen1'].chr_db == db.resolve('gen1/chr')
      assert result['gen2'].chr_db == db.resolve('gen2/chr')
      assert result['gen3'].chr_db == db.resolve('gen3/chr')

      assert result['gen1'].blast_db == db.resolve('gen1/x-blast-db')
      assert result['gen2'].blast_db == db.resolve('gen2/x-blast-db')
      assert result['gen3'].blast_db == db.resolve('gen3/x-blast-db')

}

def void testParseGenomesFolder() {

  def root = file('test-folder')

  try {
      // create the structure to test
      def folder1 = root.resolve( 'alpha')
      def folder2 = root.resolve( 'beta')
      def folder3 = root.resolve( 'delta')
      folder1.mkdirs()
      folder2.mkdirs()
      folder3.mkdirs()

      folder1.resolve('gen1.fa').text = 'uno'
      folder2.resolve('gen2.fa').text = 'due'
      folder3.resolve('gen3.fa').text = 'tre'

      def db = ('db')

      // call the function to test
      def result = parseGenomesFolder(db, root, 'y-blast')

      // verify result
      assert result.size() == 3

      assert result['alpha'].genome_fa == folder1.resolve('gen1.fa')
      assert result['beta'].genome_fa == folder2.resolve('gen2.fa')
      assert result['delta'].genome_fa == folder3.resolve('gen3.fa')

      assert result['alpha'].chr_db == db.resolve('alpha/chr')
      assert result['beta'].chr_db == db.resolve('beta/chr')
      assert result['delta'].chr_db == db.resolve('delta/chr')

      assert result['alpha'].blast_db == db.resolve('alpha/y-blast-db')
      assert result['beta'].blast_db == db.resolve('beta/y-blast-db')
      assert result['delta'].blast_db == db.resolve('delta/y-blast-db')

  }
  finally {
     root.deleteDir()
  }


}
