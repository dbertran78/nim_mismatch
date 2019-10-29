import hts
import hts/private/hts_concat

import parseopt

proc get_number_of_match(cig:Cigar): int =
  #echo cig
  for operation in cig:
    if op(operation) == match :
      return len(operation)


proc get_avg_mismatch_ratio(abam:Bam) =
  #Get the number of mismatches per read
  var nb_match, nb_indel, mismatch_total, match_total, mismatch_in_alignment, nb_read = 0
  var cig:Cigar
  for record in abam:
    #from_string*(r:Record, record_string:string) =

    #Skip unwanted reads
    if not proper_pair(flag(record)): continue
  
    #echo record,  "|", flag(record), "|"
    match_total = match_total + get_number_of_match(record.cigar)

    nb_match =0
    nb_indel = 0
    cig = cigar(record)
    for operation in cig:
      #echo op(operation)
      if op(operation) == match :
        nb_match = nb_match + 1
      if op(operation) == CigarOp.insert :
        nb_indel = nb_indel + len(operation)
      if op(operation) == CigarOp.deletion :
        nb_indel = nb_indel + len(operation)
         
    var mismatches = tag[int](record, "NM")#record.aux("NM")
    mismatch_in_alignment = mismatches.get - nb_indel#rg.tostring()
    mismatch_total = mismatch_total + mismatch_in_alignment 
    nb_read = nb_read + 1
  
  echo "Average mismatch ratio: " , "Match:", match_total , "\tMismatch:" , mismatch_total , "\tRatio:" , (mismatch_total / match_total)



### MAIN COMMANDS
var bam_file:string
var ref_file:string

#Get the arguments  
var args = initOptParser()#cmdLine)
for kind, key, val in args.getopt():
  case kind
  of cmdEnd: doAssert(false)  # Doesn't happen with getopt()
  of cmdShortOption, cmdLongOption:
    case key 
        of "bam":
          echo "Load_bam_file\t", val
          bam_file = val
        of "ref":
          echo "Load_reference_file ", val
          ref_file = val
        else:
          echo "Option unknown ", val, " ", key
          #****ADD exit code
          
  of cmdArgument:
    echo "Argument: ", key

var abam:Bam
open(abam, bam_file, index=true, fai=ref_file)


#Get the mismatc ratio
get_avg_mismatch_ratio(abam)
