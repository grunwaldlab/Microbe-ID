function validateFasta(its_blast)
{
	/*jshint multistr:true */
	var elementsInputs;
	elementsInputs = its_blast.fasta_file.value;
	var j=0;
	for(var i = 0; i < elementsInputs.length; i++) {
		if (elementsInputs[i] == '>') {
			j++;
			}
	}
	if (elementsInputs[1] == "" ){
		alert ('No sequences aaaa')
		return false;	
		}
	if (elementsInputs == "" ){
		alert ('No sequence found')
		return false;	
		}
	else if (j > 1){
		alert('Only one FASTA sequence is permitted');
		return false;
	} else if (j == 0){	
		alert('Only FASTA format is permitted');
		return false;
	}
}



function validateFasta(cox_blast)
{
	/*jshint multistr:true */
	var elementsInputs;
	elementsInputs = cox_blast.fasta_file.value;
	var j=0;
	for(var i = 0; i < elementsInputs.length; i++) {
		if (elementsInputs[i] == '>') {
			j++;
		}
	}
	
	if (elementsInputs == "" ){
			alert ('No sequence found')
			return false;	
	}	
	else if (j > 1){
		alert('Only one FASTA sequence is permitted');
		return false;
	} 
	else if (j == 0){	
		alert('Only FASTA format is permitted');
		return false;
	}	
}
