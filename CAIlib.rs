use std::process::Command;
use std::process;
use std::fs::File;
use std::str::FromStr;
use std::error::Error;
use std::{str,env};
use regex::{Regex, RegexBuilder};
use math::mean;
use std::vec::Vec;
use std::collections::HashMap;
use std::iter::Enumerate;
use std::io::{self, prelude::*, BufReader};
use std::collections::BTreeMap;
use core::ops::Range;
use std::convert::AsRef;
use std::path::Path;
use core::cmp::{min, max};


fn generate_CAI(sequences: Vec<String>, synonymous_codons: &HashMap<&str,Vec<&str>>) -> HashMap<String,f32> {
   //create a hashmap of expected sequence weights per codon as float
   let mut sequence_weights: HashMap<String, f32> = HashMap::new();
   let mut seqvec: Vec<String> = Vec::new();
   let mut seqylen: f32 = 0.0;
   for seq in sequences {
      let mut seq = &seq.to_ascii_uppercase();
      for i in 0..seq.len() {
         if i < seq.len()-3 {
            seqvec.push(seq[i..i+3].to_string());
	    }
	 }
      }
    let sv = seqvec.clone();
    for j in 0..seqvec.len() {
          *sequence_weights.entry(sv[j].to_string()).or_insert(0.0)+=1.0;
	  seqylen+=1.0;
	  }
    let numcodons = sequence_weights.keys().len();
    let seqcodons: Vec<_> = sequence_weights.keys().cloned().collect();
    for seqcodon in seqcodons {
        let val = sequence_weights[&seqcodon] as f32/(seqylen as f32/numcodons as f32);
	let val2 = synonymous_codons[&*seqcodon].len() as f32;
	let vall = if val2 != 0.0 { 1.0/val2 } else { 0.5 };
        sequence_weights.insert(seqcodon,val/vall);
	}
   sequence_weights
}

fn generate_GC(seq: String, gccontent: f32, synonymous_codons: &HashMap<&str,Vec<&str>>) -> f32 {
   //finds the GC content of a coding frame sequence as a string and returns a float of the percentage GC content
   let mut resultGC: f32 = 0.0;
   let frames = vec![0,1,2];
   for frame in frames {
       let mut cods: Vec<_> = Vec::new();
       for i in frame..seq.len() {
           if i < seq.len() - (3+frame) {
	      cods.push(seq[i..i+3].to_string());
	      }
	   }
       let mut countgc = 0;
       for cod in &cods {
           match cod.chars().last().unwrap() {
	   'G' => countgc+=1,
	   'C' => countgc+=1,
	   _   => (),
	   };
       }
       resultGC = (countgc as f32)/gccontent;
       //resultGC.insert(frame as u32,(countgc as f32)/gccontent);
       }
   resultGC
}

fn load_CAI(file: String) -> HashMap<String, f32> {
   //for loading a file with CAI stats
   let filer = File::open(&file).expect("problem opening CAI stats filename");
   let mut loaded_weights: HashMap<String,f32> = HashMap::new();
   let mut readery = BufReader::new(filer);
   for line in readery.lines() {
       let lin = line.unwrap();
       if lin.len() > 1 {
           let mut elements = lin.split(",").collect::<Vec<&str>>();
           loaded_weights.insert(elements[0].to_string(),elements[1].parse::<f32>().unwrap());
	   }
   }
   loaded_weights
}
   

fn relative_adaptiveness(seq: String, sequence_weights: HashMap<String,f32>, synonymous_codons: &HashMap<&str,Vec<&str>> ) -> HashMap<String,f32> {
   // the relative adaptiveness is the frequency of that codon compared to the freq of the optimal
   //create a hashmap with the relative adaptiveness values per codon on the sequence compared to the optimal
   let mut seqy = &seq.to_ascii_uppercase();
   let mut seqyvec: Vec<String> = Vec::new();
   let mut comparison_weights: HashMap<String, f32> = HashMap::new();
   let mut final_weights: HashMap<String, f32> = HashMap::new();
   let mut seqylen: f32 = 0.0;
   for i in 0..seqy.len() {
       if i < seqy.len()-3 {
           seqyvec.push(seqy[i..i+3].to_string());
	   }
       }
   for k in 0..seqyvec.len() {
       *comparison_weights.entry(seqyvec[k].to_string()).or_insert(0.0)+=1.0;
       seqylen+=1.0;
       }
   let seqycodons: Vec<_> = comparison_weights.keys().cloned().collect();
   let numcodons = sequence_weights.keys().len();
   for seqycodon in seqycodons {
      let valy = comparison_weights[&seqycodon] as f32/(seqylen as f32/numcodons as f32);
      let val2 = synonymous_codons[&*seqycodon].len() as f32;
      let vall = if val2 != 0.0 { 1.0/val2 } else { 0.5 };
      comparison_weights.insert(seqycodon,valy/vall);
      }
   let realcodons: Vec<_> = comparison_weights.keys().cloned().collect();
   for realcod in realcodons {
      let orig_weight = sequence_weights[&realcod];
      let new_weight = comparison_weights[&realcod];
      final_weights.insert(realcod,new_weight/orig_weight);
      }
   final_weights
}
   
   
fn CAI(testseq: String, final_weights: HashMap<String,f32>) -> f64 {
   let seqycodons: Vec<_> = final_weights.keys().cloned().collect();
   let mut collvec = Vec::new();
   //collect weights on codon sequences excluding stops and starts
   for seq in seqycodons {
      if &seq != "ATG|TGG|TGA|TAG|TAA" {
          collvec.push(final_weights[&seq] as f64);
          }
      }
   //codon adaptation index is the geometric mean of the values
   let res = mean::geometric(&collvec.as_slice());
   res
}

pub fn achieve_CAI(testing_cai: HashMap<(i32,i32),String>, gccontent: f32) -> HashMap<(i32,i32),f32> {
   //returns a hashmap of testing sequences by location and single float value of the CAI
   let mut seqsannot: AnnotMap<String, String> = AnnotMap::new();
   let mut synonymous_codons = HashMap::new();
   //inserting the codon synonymns
   synonymous_codons.insert("TTA",vec!["TTG","CTT","CTC","CTA","CTG"]);
   synonymous_codons.insert("ATG",vec!["ATG"]);
   synonymous_codons.insert("TTT",vec!["TTC"]);
   synonymous_codons.insert("TTC",vec!["TTT"]);
   synonymous_codons.insert("TAA",vec!["TGA","TAG"]);
   synonymous_codons.insert("TGA",vec!["TAA","TAG"]);
   synonymous_codons.insert("TAG",vec!["TAA","TGA"]);
   synonymous_codons.insert("TTG",vec!["TTA","CTT","CTC","CTA","CTG"]);
   synonymous_codons.insert("CTT",vec!["TTA","TTG","CTC","CTA","CTG"]);
   synonymous_codons.insert("CTC",vec!["TTA","TTG","CTT","CTA","CTG"]);
   synonymous_codons.insert("CTA",vec!["TTA","TTG","CTT","CTC","CTG"]);
   synonymous_codons.insert("CTG",vec!["TTA","TTG","CTT","CTC","CTA"]);
   synonymous_codons.insert("ATT",vec!["ATC","ATA"]);
   synonymous_codons.insert("ATC",vec!["ATT","ATA"]);
   synonymous_codons.insert("ATA",vec!["ATC","ATT"]);
   synonymous_codons.insert("GTT",vec!["GTC","GTA","GTG"]);
   synonymous_codons.insert("GTA",vec!["GTC","GTT","GTG"]);
   synonymous_codons.insert("GTC",vec!["GTT","GTA","GTG"]);
   synonymous_codons.insert("GTG",vec!["GTT","GTC","GTA"]);
   synonymous_codons.insert("TCT",vec!["TCC","TCA","TCG","AGT","AGC"]);
   synonymous_codons.insert("TCC",vec!["TCT","TCA","TCG","AGT","AGC"]);
   synonymous_codons.insert("TCA",vec!["TCT","TCC","TCG","AGT","AGC"]);
   synonymous_codons.insert("TCG",vec!["TCC","TCT","TCA","AGT","AGC"]);
   synonymous_codons.insert("AGT",vec!["TCC","TCT","TCA","TCG","AGC"]);
   synonymous_codons.insert("AGC",vec!["TCC","TCT","TCA","TCG","AGT"]);
   synonymous_codons.insert("CCT",vec!["CCC","CCA","CCG"]);
   synonymous_codons.insert("CCC",vec!["CCT","CCA","CCG"]);
   synonymous_codons.insert("CCA",vec!["CCC","CCT","CCG"]);
   synonymous_codons.insert("CCG",vec!["CCC","CCT","CCA"]);
   synonymous_codons.insert("ACT",vec!["ACC","ACA","ACG"]);
   synonymous_codons.insert("ACC",vec!["ACT","ACA","ACG"]);
   synonymous_codons.insert("ACA",vec!["ACT","ACC","ACG"]);
   synonymous_codons.insert("ACG",vec!["ACT","ACC","ACA"]);
   synonymous_codons.insert("GCT",vec!["GCC","GCA","GCG"]);
   synonymous_codons.insert("GCC",vec!["GCT","GCA","GCG"]);
   synonymous_codons.insert("GCA",vec!["GCT","GCC","GCG"]);
   synonymous_codons.insert("GCG",vec!["GCT","GCC","GCA"]);
   synonymous_codons.insert("TAT",vec!["TAC"]);
   synonymous_codons.insert("TAC",vec!["TAT"]);
   synonymous_codons.insert("CAT",vec!["CAC"]);
   synonymous_codons.insert("CAC",vec!["CAT"]);
   synonymous_codons.insert("CAA",vec!["CAG"]);
   synonymous_codons.insert("CAG",vec!["CAA"]);
   synonymous_codons.insert("AAT",vec!["AAC"]);
   synonymous_codons.insert("AAC",vec!["AAT"]);
   synonymous_codons.insert("AAA",vec!["AAG"]);
   synonymous_codons.insert("AAG",vec!["AAA"]);
   synonymous_codons.insert("GAT",vec!["GAC"]);
   synonymous_codons.insert("GAC",vec!["GAT"]);
   synonymous_codons.insert("GAA",vec!["GAG"]);
   synonymous_codons.insert("GAG",vec!["GAA"]);
   synonymous_codons.insert("TGT",vec!["TGC"]);
   synonymous_codons.insert("TGC",vec!["TGT"]);
   synonymous_codons.insert("TGG",vec!["TGG"]);
   synonymous_codons.insert("CGT",vec!["CGC","CGA","CGG","AGA","AGG"]);
   synonymous_codons.insert("CGC",vec!["CGT","CGA","CGG","AGA","AGG"]);
   synonymous_codons.insert("CGA",vec!["CGT","CGC","CGG","AGA","AGG"]);
   synonymous_codons.insert("CGG",vec!["CGT","CGC","CGA","AGA","AGG"]);
   synonymous_codons.insert("AGA",vec!["CGT","CGC","CGA","CGG","AGG"]);
   synonymous_codons.insert("AGG",vec!["CGT","CGC","CGA","CGG","AGA"]);
   synonymous_codons.insert("AGT",vec!["AGC"]);
   synonymous_codons.insert("AGC",vec!["AGT"]);
   synonymous_codons.insert("GGT",vec!["GGC","GGA","GGG"]);
   synonymous_codons.insert("GGC",vec!["GGT","GGA","GGG"]);
   synonymous_codons.insert("GGA",vec!["GGT","GGC","GGG"]);
   synonymous_codons.insert("GGG",vec!["GGT","GGA","GGC"]);
   let mut list_o_seqs: Vec<String> = Vec::new();
   for (_,seqa) in &testing_cai {
      list_o_seqs.push(seqa.to_string());
      }
   let seqkeys: Vec<_> = testing_cai.keys().cloned().collect();
   let mut final_result: f32 = 0.0;
   let mut return_results: HashMap<(i32,i32),f32> = HashMap::new();
   //collect the relative adaptiveness values per tested sequence
   let mut results = generate_CAI(list_o_seqs, &synonymous_codons);
   for seqkey in seqkeys {
     let newi = relative_adaptiveness(testing_cai[&seqkey].clone(),results.clone(),&synonymous_codons);
     final_result = CAI(testing_cai[&seqkey].clone(), newi);
     return_results.insert(seqkey,final_result);
     }
   return_results
}
