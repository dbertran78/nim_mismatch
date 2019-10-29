import hts
#import hts/private/hts_concat
import strutils
import parseopt

type
  StringArray = array[0..4, string] # an array that is indexed with 0..5
const num_base_pair: StringArray = ["A", "T", "G", "C", "N"]
type
  mat5x5 = array[0..4, array[0..4, int]]
  
proc convert_base_pair(bp:char): int =
  case bp:
    of 'A':
      return 0
    of 'T':
      return 1
    of 'G':
      return 2
    of 'C':
       return 3
    of 'N':
      return 4
    else:
      echo "Unknow base pair bp"
      return -1
      
  
proc update_bam_cigar(
  abam:Bam,
  obam:var Bam,
  fai:Fai,
  mut_type:var mat5x5
     ) =
  var ref_seq:string
  var chr:string
  var r_start, r_end = 0
  var nb_mismatch,read_counter, ref_counter, equal_counter, x_counter = 0
  var read_seq:string
  var cig:Cigar
  var new_cigar:string
  var new_string_record: seq[string]
  for record in abam:
    
    #from_string*(r:Record, record_string:string) =
    #var s:string = "5821:596	99	1000000	596	60	100M	=	896	400	TTTTTCCCTTAATCCCAATCATTATGAGAAGTCTGTAGGTCAAGTGATACACAAATAGCACAGGACCGGTAAGCTAGTACGCATCTGATTTCTGAGCCGG	*	NM:i:0	MD:Z:100	MC:Z:100M	AS:i:100	XS:i:0"
    var aln_info = tag[string](record, "MD")#record.aux("NM")

    #Skip unwanted reads
    if not proper_pair(flag(record)): continue

    chr = chrom(record)
    r_start = start(record)
    r_end = stop(record)
    cig = cigar(record)
    #echo qname(record), " ", aln_info, " ", cig, " ", chr, " ", r_start, " ", r_end
    
    ref_seq=get(fai, chrom(record), r_start, r_end-1)
    sequence(record, read_seq)
    nb_mismatch = 0
    read_counter=0
    ref_counter = 0
    x_counter = 0
    equal_counter = 0
    new_cigar = ""
    
    #echo ref_seq , "\n", read_seq, "\n", new_cigar
    #Loop threw the cigar string to identify the mismatches
    for operation in cig:
      #Match
      if op(operation) == match :
        #Loop threw ref and read and find mismatches
        for i in 0..len(operation)-1:
          #This is a mismatch
          if ref_seq[ref_counter + i] != read_seq[read_counter + i] and ref_seq[ref_counter + i] != 'N' :
            #echo "mismatch  ", (ref_counter + i), " ", ref_seq[ref_counter + i], " ", read_seq[read_counter + i]
            #echo convert_base_pair(ref_seq[ref_counter + i]) , " ", convert_base_pair(read_seq[read_counter + i])
            mut_type[convert_base_pair(ref_seq[ref_counter + i])][convert_base_pair(read_seq[read_counter + i])] = mut_type[convert_base_pair(ref_seq[ref_counter + i])][convert_base_pair(read_seq[read_counter + i])] + 1
            nb_mismatch = nb_mismatch+1
            x_counter = x_counter + 1
            if equal_counter != 0:
              new_cigar = new_cigar & $equal_counter & "="
              equal_counter = 0
          else:
            equal_counter = equal_counter + 1
            if x_counter != 0:
              new_cigar = new_cigar & $x_counter & "X"
              x_counter = 0
        #Print the last operation in the new cigar
        if x_counter != 0:
          new_cigar = new_cigar & $x_counter & "X"
        else:
          new_cigar = new_cigar & $equal_counter & "="
      else:
        new_cigar = new_cigar & $operation
        equal_counter=0
        x_counter=0

      #Update the counter base on the non-match cigar string operation
      if query(consumes(operation)):
        read_counter = read_counter + len(operation)

      if reference(consumes(operation)):
        ref_counter = ref_counter + len(operation)

    #var mismatches = tag[int](record, "NM")#record.aux("NM")
    #var mismatch_in_alignment = mismatches.get#rg.tostring()

    #Create the new bam file
#from_string(record, s)
#    obam.write(record)
#    obam.write(record)
        
    ######################################3
    var nb_indel = 0
    for operation in cig:
      #echo op(operation)
      if op(operation) == CigarOp.insert :
        nb_indel = nb_indel + len(operation)
      if op(operation) == CigarOp.deletion :
        nb_indel = nb_indel + len(operation)
        
    var mismatches = tag[int](record, "NM")#record.aux("NM")
    var mismatch_in_alignment = mismatches.get - nb_indel
    if nb_mismatch != mismatch_in_alignment :
      echo qname(record), " ", aln_info, " ", cig, " ", chr, " ", r_start, " ", r_end, " ",  new_cigar, " ", read_counter, " ",  ref_counter
      echo ref_seq , "\n", read_seq, "\n", nb_mismatch , "!=",  mismatch_in_alignment, " ", nb_indel

    ##################################
      
    #echo ref_seq , "\n", read_seq, "\n", new_cigar
    new_string_record = split(to_string(record))#new_string_record =
    #echo "\n", new_string_record
    new_string_record[5] = new_cigar
    #echo new_string_record
    #  qname(record) , "\t",
      
    from_string(record, new_string_record.join("\t"))
    obam.write(record)
    
  
### MAIN COMMANDS
var bam_file:string
var sam_file:string
var ref_file:string
var mut_type:mat5x5
for i in 0..4:
  for j in 0..4:
    mut_type[i][j] = 0

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
        of "new-sam":
          echo "Write_sam_file\t", val
          sam_file = val
        of "ref":
          echo "Load_reference_file ", val
          ref_file = val
        else:
          echo "Option unknown ", val, " ", key
          #****ADD exit code
          
  of cmdArgument:
    echo "Argument: ", key

var abam:Bam
open(abam, bam_file, index=true)#, fai=ref_file)
var fai:Fai
if not open(fai, ref_file):
  echo "Problem opening ", ref_file

#
var obam:Bam
open(obam, sam_file, mode="w")
obam.write_header(abam.hdr)
update_bam_cigar(abam, obam, fai, mut_type)

echo "\nmismatch_types: "
var nb_mismatch = 0
for i in 0..4:
  for j in 0..4:
    if i != j and mut_type[i][j] != 0:
      nb_mismatch = nb_mismatch + mut_type[i][j]
      echo num_base_pair[i], "->", num_base_pair[j], " ", mut_type[i][j]

echo "\nnb_mismatch: ", nb_mismatch, "\n"

obam.close()




