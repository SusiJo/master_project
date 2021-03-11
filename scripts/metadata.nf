/*
* usage:
* nextflow run metadata.nf --run_acc_list <path_to_run_accs> --outdir results
*/

params.experiment = 'https://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?save=efetch&db=sra&rettype=SraExperimentPackage&term='
params.runinfo = 'https://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?save=efetch&db=sra&rettype=runinfo&term='

/*
 * Step 1 - prefetch
 */ 

process prefetch {
    input:
    val run_acc from Channel.fromPath(params.run_acc_list).splitText()

    output:
    val output_file into ch_accessions // needs to be added to prefetch process

    script:
    output_file = run_acc.trim()
    """
    """    
}

/*
 * Step 5 - Download Metadata XML, CSV
 */

process get_metadata {
	publishDir "${params.outdir}/metadata", mode:'copy'  	

	input:
	val run_acc from ch_accessions 		

	output:
	//file "[S,E,D]RR*[0-9].xml"
	file "[S,E,D]RR*[0-9].csv"		

	script:
        def acc = run_acc
	//def experiment_url = "${params.experiment}" + "${acc}"	wget '${experiment_url}' -nv -O ${acc}.xml curl '${experiment_url}' -s --output ${acc}.xml
	def runinfo_url = "${params.runinfo}" + "${acc}"	
	"""	
	if [ -x /usr/bin/wget ] ; then		
		wget '${runinfo_url}' -nv -O ${acc}.csv
	else
		curl '${runinfo_url}' -s --output ${acc}.csv
	fi	
	"""	
}
