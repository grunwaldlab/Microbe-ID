function validateForm()
{
	/*jshint multistr:true */
var elementsInputs=document.forms["fasta_file"]["fasta"].value;
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
	else if (j == 0){	
		alert('Only FASTA format is permitted');
		return false;
	}

}