function validateForm()
{
var x=document.fasta_file.fasta.value;
var atpos=x.indexOf(">");
if (atpos<1)
  {
  alert("Not a valid FASTA file");
  return false;
  }
}