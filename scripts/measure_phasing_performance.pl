#!/usr/bin/env perl

use strict;
use warnings;

use vars qw/ %opt /;
use Class::Struct;
use IO::Zlib;

# configs

sub init()
{
	use Getopt::Std;
	my $opt_string = 'hi:r:o:e:a:';
	getopts( "$opt_string", \%opt ) or usage();
	usage() if ($opt{h} || !$opt{i} || !$opt{r} || !$opt{o});
}

sub usage()
{
	print STDERR "usage: $0 [-h] -i input_vcf_file_to_be_evaluated -r reference_file(.vcf, SNP only) -o output_file_prefix [-a output_for_pairwise_SNV_analysis]\n";
	exit;
}

init();


struct SNP =>
[
	# golden set
	_pos => '$',	# position
	_phased => '$',	# 1: phased, 0: not phased
	_nt1 => '$',	# first nucleotide
	_nt2 => '$',	# second nucleotide
	_type => '$',	# "ref", "het", or "hom"
	_nt_ref => '$',	# reference nucleotide
	
	# phasing prediction
	_r_block_id => '*$',	# 0: no info (unphased), other integers: block id
	_r_phase_info => '*$',	# 1: same as golden, -1: opposite, 0: no information or incorrect allele
	_r_var_info => '*$',	# 0: no variant call, -1: not matched, 1: matched

	# for performance measure
	_r_first_in_block => '*$',	# 1: yes
	_r_error => '*$',			# 0:correct, 1:undetermined_switch_error, 2:long_switch_error, 3:point_switch_error, 4:first_pos, 5:unphased, 6:tmp_switch_error
	_r_last_in_block => '*$',	# 1: yes
];

my %count = ();
my %count_desc = ();
my %snps = ();
my $chr_check = "";


my $desc = << 'DESCRIPTION';
blocks.measured	# Phasing blocks considered for evaluation
blocks.unmeasured	# Phasing blocks not considered for evaluation (because of no phased sites)
num.blocks	# Phasing blocks in the VCF evaluated (simply by PS tags)
pred.het	# Heterozygous sites in the VCF evaluated
pred.het.phased	# Phased heterozygous sites in the VCF evaluated
pred.hom	# Homozygous sites in the VCF evaluated
pred.not_in_ref	# Sites in the VCF evaluated, not represented in the golden set
pred.not_in_ref.het	# Heterozygous sites in the VCF evaluated, not represented in the golden set
pred.not_in_ref.het.phased	# Phased heterozygous sites in the VCF evaluated, not represented in the golden set
pred.not_in_ref.hom	# Homozygous sites in the VCF evaluated, not represented in the golden set
pred.total	# Total sites in the VCF evaluated
ref.het	# Heterozygous sites in the golden set
ref.het.phased	# Phased heterozygous sites in the golden set
ref.het.phased.correct_call_in_pred	# Phased heterozygous sites in the golden set that have consistent genotypes calls in the VCF evaluated
ref.het.phased.correct_call_in_pred.phased_in_pred	# Phased heterozygous sites in the golden set that have consistent phased genotypes in the VCF evaluated
ref.het.phased.correct_call_in_pred.phased_in_pred.long_switch_error	# Long switch errors
ref.het.phased.correct_call_in_pred.phased_in_pred.point_switch_error	# Point switch errors
ref.het.phased.correct_call_in_pred.phased_in_pred.switch_error	# Total switch errors (long + 2*switch + undetermined)
ref.het.phased.correct_call_in_pred.phased_in_pred.switch_error_measured	# Phased heterozygous sites in the golden set that used for the evaluation
ref.het.phased.correct_call_in_pred.phased_in_pred.undetermined_switch_error	# Switch errors not determined (not long nor switch)
ref.het.phased.incorrect_call_in_pred	# Phased heterozygous sites in the golden set where genotypes in the VCF evaluated are different
ref.het.phased.no_call_in_pred	# Phased heterozygous sites in the golden set where no genotypes calls were made in the VCF evaluated
ref.het.unphased	# Unphased heterozygous sites in the golden set
ref.hom	# Homozygous sites in the golden set
ref.total	# Total sites in the golden set
DESCRIPTION


#	Open std and error output files (storing switch error positions)
open(my $fh_o, ">$opt{o}.out") or die "cannot open file: $!";
open(my $fh_e, ">$opt{o}.err_pos") or die "cannot open file: $!";
open(my $fh_e2, ">$opt{o}.err_pos.more") or die "cannot open file: $!";

foreach my $line (split/\n/, $desc) {
	chomp $line;
	my ($key, $d) = split/\t/, $line;
	set_counter_desc($key, $d);
}

# Load golden standard set
load_golden_set($opt{r});
# Load phased result from a phasing method
my $phasing_block_num = read_phasing_block_vcf($opt{i});


####
# Statistics 
####
my @pos_sorted = sort {$a<=>$b} keys %snps;
my @pos_phased_sorted;
foreach my $k (@pos_sorted) {
	if ($snps{$k}->_type eq "het" && $snps{$k}->_phased == 1) {
		push(@pos_phased_sorted, $k);
	}
}

####
# Basic counts
####
foreach my $k (@pos_sorted) {
	inc_count("ref.total");
	if ($snps{$k}->_type eq "het") {
		inc_count("ref.het");
		if ($snps{$k}->_phased == 1) {
			inc_count("ref.het.phased");
			if (${$snps{$k}->_r_var_info} == 0) {
				inc_count("ref.het.phased.no_call_in_pred");
			} elsif (${$snps{$k}->_r_var_info} == -1) {
				inc_count("ref.het.phased.incorrect_call_in_pred");
			} else {
				inc_count("ref.het.phased.correct_call_in_pred");
				if (${$snps{$k}->_r_phase_info} != 0) {
					inc_count("ref.het.phased.correct_call_in_pred.phased_in_pred");
				} else {
					inc_count("ref.het.phased.correct_call_in_pred.unphased_in_pred");
				}
			}
		} else {
			inc_count("ref.het.unphased");
		}
	} elsif ($snps{$k}->_type eq "hom") {
		inc_count("ref.hom");
	} elsif ($snps{$k}->_type eq "ref") {
		inc_count("ref.ref");
	} else {
		inc_count("ref.unknown_type");
	}
}

printf $fh_o "- Reference set statistics\n";
printf $fh_o "# Total SNPs in the reference set:\t%d\n", scalar keys %snps;
printf $fh_o "# Phasing blocks in the prediction:\t%d\n", $phasing_block_num;

####
# SER
####
my $cutoff_dist_for_long_switch = 3;
my @last_phase = (0)x($phasing_block_num+1);
my @last_phase_long = (0)x($phasing_block_num+1);
my @dist_from_last_switch_error = (0)x($phasing_block_num+1);	# in number of phased variants
my @dist_from_first_pos = (0)x($phasing_block_num+1);	# in number of phased variants
my @last_phased_pos = (-1)x($phasing_block_num+1);
my @last_undetermined_switch_error_pos = (-1)x($phasing_block_num+1);
foreach my $k (@pos_phased_sorted) {
	my $b = ${$snps{$k}->_r_block_id};
	if ($b == 0) {
		# no information from input file
		${$snps{$k}->_r_error} = 5;
		next;
	}
	if (${$snps{$k}->_r_phase_info} != 0) {
		if ($last_phase[$b] != 0) {
			$dist_from_first_pos[$b]++;
			inc_count("ref.het.phased.correct_call_in_pred.phased_in_pred.switch_error_measured");
			if ($last_phase[$b] != ${$snps{$k}->_r_phase_info}) {
				inc_count("ref.het.phased.correct_call_in_pred.phased_in_pred.switch_error");
				my $dist_to_last_phased_pos = $k - $last_phased_pos[$b];
				printf $fh_e "%d\t%d\n", $snps{$k}->_pos, $dist_to_last_phased_pos;

				if ($last_undetermined_switch_error_pos[$b] > 0) {
					if ($dist_from_last_switch_error[$b] < $cutoff_dist_for_long_switch) {
						# considered as a point switch
						${$snps{$last_undetermined_switch_error_pos[$b]}->_r_error} = 3;
						${$snps{$k}->_r_error} = 3;
						$last_undetermined_switch_error_pos[$b] = -1;
						inc_count("ref.het.phased.correct_call_in_pred.phased_in_pred.point_switch_error");
						dec_count("ref.het.phased.correct_call_in_pred.phased_in_pred.undetermined_switch_error");
					} else {
						# considered as an undetermined switch for now
						${$snps{$k}->_r_error} = 6;
						$last_undetermined_switch_error_pos[$b] = $k;
						inc_count("ref.het.phased.correct_call_in_pred.phased_in_pred.undetermined_switch_error");
					}
				} else {
					# considered as an undetermined switch for now
					if ($dist_from_first_pos[$b] < $cutoff_dist_for_long_switch) {
						${$snps{$k}->_r_error} = 1;
					} else {
						${$snps{$k}->_r_error} = 6;
					}
					$last_undetermined_switch_error_pos[$b] = $k;
					inc_count("ref.het.phased.correct_call_in_pred.phased_in_pred.undetermined_switch_error");
				}
				$dist_from_last_switch_error[$b] = 0;
			} else {
				${$snps{$k}->_r_error} = 0;
				$dist_from_last_switch_error[$b]++;
				if ($last_undetermined_switch_error_pos[$b] > 0 && $dist_from_last_switch_error[$b] == $cutoff_dist_for_long_switch) {
					if (${$snps{$last_undetermined_switch_error_pos[$b]}->_r_error} == 6) {
						# last error is considered as a long switch
						${$snps{$last_undetermined_switch_error_pos[$b]}->_r_error} = 2;
						dec_count("ref.het.phased.correct_call_in_pred.phased_in_pred.undetermined_switch_error");
						inc_count("ref.het.phased.correct_call_in_pred.phased_in_pred.long_switch_error");
						$last_undetermined_switch_error_pos[$b] = -1;
					}
				}
			}
		} else {
			${$snps{$k}->_r_error} = 4;
		}	
		$last_phase[$b] = ${$snps{$k}->_r_phase_info};
		$last_phased_pos[$b] = $k;
	} else {
		${$snps{$k}->_r_error} = 5;
	}
}

foreach my $k (@pos_phased_sorted) { 
	my $msg;
	if (${$snps{$k}->_r_error} == 0) {
		$msg = "CORRECT";
	} elsif (${$snps{$k}->_r_error} == 1 || ${$snps{$k}->_r_error} == 6) {
		$msg = "SWITCH_ERROR_UNDEF";
	} elsif (${$snps{$k}->_r_error} == 2) {
		$msg = "SWITCH_ERROR_LONG";
	} elsif (${$snps{$k}->_r_error} == 3) {
		$msg = "SWITCH_ERROR_POINT";
	} elsif (${$snps{$k}->_r_error} == 4) {
		$msg = "FIRST_POS";
	} elsif (${$snps{$k}->_r_error} == 5) {
		$msg = "UNPHASED";
	} else {
		my $code = ${$snps{$k}->_r_error};
		printf STDERR "ERROR: undefined error code ($code)\n";
		exit;
	}
	printf $fh_e2 "%d\t%s\n", $snps{$k}->_pos, $msg; 
}

close($fh_e);
close($fh_e2);

printf $fh_o "SER(%%):\t%.2f\n", get_count("ref.het.phased.correct_call_in_pred.phased_in_pred.switch_error")*100/get_count("ref.het.phased.correct_call_in_pred.phased_in_pred.switch_error_measured");

####
# N50 (this is based on the subset overlapped with golden set)
####
my @first_pos = (0)x($phasing_block_num+1);	# first phased position
my @last_pos = (0)x($phasing_block_num+1);	# last phased position
foreach my $k (@pos_phased_sorted) {
	my $b = ${$snps{$k}->_r_block_id};
	if (${$snps{$k}->_r_phase_info} != 0) {
		if ($first_pos[$b] == 0) {
			$first_pos[$b] = $k;
		}
		$last_pos[$b] = $k;
	}
}
my @block_lengths = (0)x($phasing_block_num+1);
my $sum_lengths = 0;
for (my $b=1; $b<=$phasing_block_num; $b++) {
	if ($first_pos[$b] == 0) {
		inc_count("blocks.unmeasured");
		next;
	}
	inc_count("blocks.measured");
	$block_lengths[$b] = $last_pos[$b] - $first_pos[$b] +1;
	$sum_lengths += $block_lengths[$b];
}
my @sorted_block_lengths = sort {$b <=> $a} @block_lengths;
my $running_sum = 0;
my $N50 = -1;
for (my $i=0; $i<$phasing_block_num; $i++) {
	$running_sum += $sorted_block_lengths[$i];
	if ($running_sum >= $sum_lengths*0.5) {
		$N50 = $sorted_block_lengths[$i];
		last;
	}
}
printf $fh_o "N50(bp):\t%d\n", $N50;
####

####
# Marking first and last phased positions for each block
####
foreach my $pos (@first_pos[1..$#first_pos]) {
	if ($pos != 0) {
		if (${$snps{$pos}->_r_first_in_block} == 1) {
			printf STDERR "More than one blocks have same boundary positions!\n";
		}
		${$snps{$pos}->_r_first_in_block} = 1;
	}
}
foreach my $pos (@last_pos[1..$#last_pos]) {
	if ($pos != 0) {
		if (${$snps{$pos}->_r_last_in_block} == 1) {
			printf STDERR "More than one blocks have same boundary positions!\n";
		}
		${$snps{$pos}->_r_last_in_block} = 1;
	}
}
####

####
# For S50, N50, AN50 (Adjusted N50)
####
struct block_info =>
[
	_first_pos => '*$',
	_last_pos => '*$',
	_num_SNPs_total => '*$',	# excluding the first SNP
	_num_SNPs_phased => '*$',	# excluding the first SNP
	_num_SNPs_unphased => '*$',	
	_last_phased_pos => '*$',
];
my @blocks;
for (my $b=1; $b<=$phasing_block_num; $b++) {
	$blocks[$b] = block_info->new(_first_pos=>$first_pos[$b],
									_last_pos=>$last_pos[$b],
									_num_SNPs_total=>0,
									_num_SNPs_phased=>0,
									_num_SNPs_unphased=>0);
}
my %cur_blocks = ();
foreach my $pos (@pos_phased_sorted) {
	foreach my $b (keys %cur_blocks) {
		${$blocks[$b]->_num_SNPs_total} += 1;
		if (${$snps{$pos}->_r_phase_info} == 0) {
			${$blocks[$b]->_num_SNPs_unphased} += 1;
		} else {
			if (${$snps{$pos}->_r_block_id} == $b) {	# phased in the block
				${$blocks[$b]->_num_SNPs_phased} += 1;
			} else {
				${$blocks[$b]->_num_SNPs_unphased} += 1;
			}
		}
	}
	if (${$snps{$pos}->_r_first_in_block} == 1) {
		$cur_blocks{${$snps{$pos}->_r_block_id}} = 1;
	}
	if (${$snps{$pos}->_r_last_in_block} == 1) {
		delete $cur_blocks{${$snps{$pos}->_r_block_id}};
	}
}

print $fh_o join("\t", "#BL, block_length(S,N,AN)", "BLOCK_ID", "NUM_SNP_TOTAL", "NUM_SNP_PHASED", "NUM_SNP_UNPHASED", "LENGTH", "ADJUSTED_LENGTH") ."\n";
for (my $b=1; $b<=$phasing_block_num; $b++) {
	if (${$blocks[$b]->_num_SNPs_total} > 0) {
		my $len = ${$blocks[$b]->_last_pos} - ${$blocks[$b]->_first_pos} +1;	# length including bases at both ends
		printf $fh_o "BL\t%d\t%d\t%d\t%d\t%d\t%.1f\n", $b, 
								${$blocks[$b]->_num_SNPs_total}+1,
								${$blocks[$b]->_num_SNPs_phased}+1,
								${$blocks[$b]->_num_SNPs_unphased},
								$len,
								$len*${$blocks[$b]->_num_SNPs_phased}/${$blocks[$b]->_num_SNPs_total};
	}
}
####

####
# For QAN50 (Quality Adjusted N50)
####
@blocks = ();
for (my $b=1; $b<=$phasing_block_num; $b++) {
	$blocks[$b] = block_info->new(_first_pos=>$first_pos[$b], _last_pos=>$last_pos[$b], _num_SNPs_total=>0, _num_SNPs_phased=>0, _num_SNPs_unphased=>0, _last_phased_pos=>$first_pos[$b]);
}
my @additional_blocks = ();
my $num_additional_blocks = 0;
%cur_blocks = ();
foreach my $pos (@pos_phased_sorted) {
	foreach my $b (keys %cur_blocks) {
		if (${$snps{$pos}->_r_phase_info} == 0) {
			${$blocks[$b]->_num_SNPs_total} += 1;
			${$blocks[$b]->_num_SNPs_unphased} += 1;
		} else {
			if (${$snps{$pos}->_r_block_id} == $b) {	# phased in the block
				if (${$snps{$pos}->_r_error} == 1 || ${$snps{$pos}->_r_error} == 2 || ${$snps{$pos}->_r_error} == 3 || ${$snps{$pos}->_r_error} == 6) {	# error; cut the phasing block and store the information
					if (${$blocks[$b]->_num_SNPs_phased} > 0) {
						$additional_blocks[$num_additional_blocks] = block_info->new(_first_pos=>${$blocks[$b]->_first_pos}, _last_pos=>${$blocks[$b]->_last_phased_pos}, _num_SNPs_total=>${$blocks[$b]->_num_SNPs_total}, _num_SNPs_phased=>${$blocks[$b]->_num_SNPs_phased}, _num_SNPs_unphased=>${$blocks[$b]->_num_SNPs_unphased});
						$num_additional_blocks++;
					}
					# re-init. for current block
					${$blocks[$b]->_first_pos} = $pos;
					${$blocks[$b]->_num_SNPs_total} = 0;
					${$blocks[$b]->_num_SNPs_phased} = 0;
					${$blocks[$b]->_num_SNPs_unphased} = 0;
				} else {
					${$blocks[$b]->_num_SNPs_total} += 1;
					${$blocks[$b]->_num_SNPs_phased} += 1;
					${$blocks[$b]->_last_phased_pos} = $pos;
				}
			} else {
				${$blocks[$b]->_num_SNPs_total} += 1;
				${$blocks[$b]->_num_SNPs_unphased} += 1;
			}
		}
	}
	if (${$snps{$pos}->_r_first_in_block} == 1) {
		$cur_blocks{${$snps{$pos}->_r_block_id}} = 1;
	}
	if (${$snps{$pos}->_r_last_in_block} == 1) {
		delete $cur_blocks{${$snps{$pos}->_r_block_id}};
	}
}

print $fh_o join("\t", "#QBL, block_length(QAN)", "BLOCK_ID", "NUM_SNP_TOTAL", "NUM_SNP_PHASED", "NUM_SNP_UNPHASED", "LENGTH", "ADJUSTED_LENGTH") ."\n";
# additional blocks
for (my $b=0; $b<$num_additional_blocks; $b++) {
	if (${$additional_blocks[$b]->_num_SNPs_total} > 0) {
		my $len = ${$additional_blocks[$b]->_last_pos} - ${$additional_blocks[$b]->_first_pos} +1;	# length including bases at both ends
		printf $fh_o "QBL\tA%d\t%d\t%d\t%d\t%d\t%.1f\n", $b, 
								${$additional_blocks[$b]->_num_SNPs_total}+1,
								${$additional_blocks[$b]->_num_SNPs_phased}+1,
								${$additional_blocks[$b]->_num_SNPs_unphased},
								$len,
								$len*${$additional_blocks[$b]->_num_SNPs_phased}/${$additional_blocks[$b]->_num_SNPs_total};
	}
}
for (my $b=1; $b<=$phasing_block_num; $b++) {
	if (${$blocks[$b]->_num_SNPs_total} > 0) {
		my $len = ${$blocks[$b]->_last_pos} - ${$blocks[$b]->_first_pos} +1;	# length including bases at both ends
		printf $fh_o "QBL\t%d\t%d\t%d\t%d\t%d\t%.1f\n", $b, 
								${$blocks[$b]->_num_SNPs_total}+1,
								${$blocks[$b]->_num_SNPs_phased}+1,
								${$blocks[$b]->_num_SNPs_unphased},
								$len,
								$len*${$blocks[$b]->_num_SNPs_phased}/${$blocks[$b]->_num_SNPs_total};
	}
}
####


# print counts
foreach my $k (sort keys %count) {
	printf $fh_o "%s:\t%d\t%s\n", $k, $count{$k}, get_counter_desc($k);
}


if (!$opt{a}) {
	exit;
}
	
####
# Phasing yield and accuracy as a function of distance
####
my $max_dist = 10000000;
my $bin_size = 1000;
my $counter_size = $max_dist/$bin_size;
my @count_phased = (0) x $counter_size;
my @count_unphased = (0) x $counter_size;
my @count_phased_correct = (0) x $counter_size;

open(my $fh_a, ">$opt{a}") or die "cannot open file: $!";
for (my $i=0; $i<$#pos_phased_sorted; $i++) {
	for (my $j=$i+1; $j<=$#pos_phased_sorted; $j++) {
		my $pos1 = $pos_phased_sorted[$i];
		my $pos2 = $pos_phased_sorted[$j];
		my $dist = $pos2 - $pos1;
		if ($dist > $max_dist) {
			last;
		}
		my $bin_idx = int(($dist-1)/$bin_size);
		my $b1 = ${$snps{$pos1}->_r_block_id};
		my $b2 = ${$snps{$pos2}->_r_block_id};
		if ($b1 == 0 || $b2 == 0) {
			# unphased
			$count_unphased[$bin_idx]++;
		} else {
			if ($b1 == $b2) {
				# same phasing block
				my $p1 = ${$snps{$pos1}->_r_phase_info};
				my $p2 = ${$snps{$pos2}->_r_phase_info};
				if ($p1 == 0 || $p2 == 0) {
					# unphased
					$count_unphased[$bin_idx]++;
				} elsif ($p1*$p2 == 1) {
					# correct
					$count_phased[$bin_idx]++;
					$count_phased_correct[$bin_idx]++;
				} elsif ($p1*$p2 == -1) {
					# incorrect
					$count_phased[$bin_idx]++;
				} else {
					# unexpected, error
					print STDERR "Unexpected values - phase=($p1, $p2), block=($b1, $b2), positions=($pos1, $pos2)\n";
					exit;
				}
			} else {
				# different phasing block
				$count_unphased[$bin_idx]++;
			}
		}
	}
}

my $cum_count_total = 0;
my $cum_count_phased = 0;
my $cum_count_phased_correct = 0;
my ($yield, $acc, $cum_yield, $cum_acc);

printf $fh_a "DISTANCE\tYIELD\tACCURACY\tYIELD_CUM\tACCURACY_CUM\tNUM_PAIR\tNUM_PAIR_SAME_BLOCK\tNUM_PAIR_CORRECT_PHASE\n";
for (my $i=0; $i<$counter_size; $i++) {
	$cum_count_total += $count_phased[$i]+$count_unphased[$i];
	$cum_count_phased += $count_phased[$i];
	$cum_count_phased_correct += $count_phased_correct[$i];

	if ($count_phased[$i]+$count_unphased[$i] == 0) {
		$yield = "NA";
	} else {
		$yield = sprintf("%.2f", $count_phased[$i]*100/($count_phased[$i]+$count_unphased[$i]));
	}

	if ($count_phased[$i] == 0) {
		$acc = "NA";
	} else {
		$acc = sprintf("%.2f", $count_phased_correct[$i]*100/$count_phased[$i]);
	}

	if ($cum_count_total == 0) {
		$cum_yield = "NA";
	} else {
		$cum_yield = sprintf("%.2f", $cum_count_phased*100/$cum_count_total);
	}
	
	if ($cum_count_phased == 0) {
		$cum_acc = "NA";
	} else {
		$cum_acc = sprintf("%.2f", $cum_count_phased_correct*100/$cum_count_phased);
	}

	printf $fh_a "%d\t%s\t%s\t%s\t%s\t%d\t%d\t%d\n", ($i+1)*$bin_size,
											$yield,
											$acc,
											$cum_yield,
											$cum_acc,
											$count_phased[$i]+$count_unphased[$i],
											$count_phased[$i],
											$count_phased_correct[$i];
}

close($fh_a);
close($fh_o);

#### MAIN end



### functions

sub dec_count {
	my $k = shift;
	if (exists($count{$k})) {
		$count{$k}--;
	} else {
		$count{$k} = 0;
	}
}

sub inc_count {
	my $k = shift;
	if (exists($count{$k})) {
		$count{$k}++;
	} else {
		$count{$k} = 1;
	}
}

sub get_count {
	my $k = shift;
	if (exists($count{$k})) {
		return $count{$k};
	} else {
		return 0;
	}
}

sub set_counter_desc {
	my $k = shift;
	my $desc = shift;
	$count_desc{$k} = $desc;
	$count{$k} = 0;
}

sub get_counter_desc {
	my $k = shift;
	if (exists($count_desc{$k})) {
		return $count_desc{$k};
	} else {
		return "No description available";
	}
}

sub read_phasing_block_vcf {
	my $file = shift;
	my %PS2blockID = ();
	my $num_PS = 0;
	my $fh;
	if ($file=~m/\.gz$/) {
		$fh = new IO::Zlib;
		$fh->open("$file", "rb") or die "cannot open file: $!";
	} else {
		open($fh, "<$file") or die "cannot open file: $!";
	}
	while(<$fh>) {
		if (m/^#/) {
			next;
		}
		chomp;
		my ($chr, $pos, $t1, $ref, $alt, $t2, $t3, $t4, $format, $v) = split/\t/, $_;
		if ($chr_check ne $chr) {
			print STDERR "chromosome name unmatched: $file, $_\n";
			exit;
		}
		if (length($ref)>1) {
			print STDERR "Non-SNP variant exists: $file, $_\n";
			exit;
		}
		my @alts = split/,/, $alt;
		for (my $i=0; $i<=$#alts; $i++) {
			if (length($alts[$i])>1) {
				print STDERR "Non-SNP variant exists: $file, $_\n";
				exit;
			}
		}
	
		my $phased = 0;
		my ($nt1, $nt2, $type);

		my @IDs = split/:/, $format;
		my @values = split/:/, $v;

		my $gt = "";
		my $ps = "";

		for (my $i=0; $i<=$#IDs; $i++) {
			if ($IDs[$i] eq "GT") {
				$gt = $values[$i];
			} elsif ($IDs[$i] eq "PS") {
				$ps = $values[$i];
				if (!exists($PS2blockID{$ps})) {
					$num_PS++;
					inc_count("num.blocks");
					$PS2blockID{$ps} = $num_PS;
				}
			}
		}
	
		if ($gt =~ m/^(\d+)([\|\/])(\d+)/) {
			my $idx1 = $1;
			my $idx2 = $3;
			if ($2 eq "|") {
				$phased = 1;
			}
			if ($idx1 == $idx2) {
				# hom
				if ($idx1 == 0) {
					$nt1 = $nt2 = $ref;
					$type = "ref";
				} else {
					$nt1 = $nt2 = $alts[$idx1-1];
					$type = "hom";
				}
			} else {
				# het
				$type = "het";
				if ($idx1 == 0) {
					$nt1 = $ref;
				} else {
					$nt1 = $alts[$idx1-1];
				}
				if ($idx2 == 0) {
					$nt2 = $ref;
				} else {
					$nt2 = $alts[$idx2-1];
				}
			}

			if ($phased == 1) {
				if ($ps eq "") {
					$ps = 0;
				}
				if (!exists($PS2blockID{$ps})) {
					$num_PS++;
					inc_count("num.blocks");
					$PS2blockID{$ps} = $num_PS;
				}
			}
		} else {
			print STDERR "unexpected format: $_\n";
			exit;
		}

		if (exists($snps{$pos})) {
			my $snp = $snps{$pos};
			if ($type eq "het") {
				if ($nt1 eq $snp->_nt1 && $nt2 eq $snp->_nt2) {
					${$snps{$pos}->_r_var_info} = 1;
					if ($phased == 1) {
						${$snps{$pos}->_r_phase_info} = 1;
					}
					if ($ps ne "" && exists($PS2blockID{$ps})) {
						${$snps{$pos}->_r_block_id} = $PS2blockID{$ps};
					}
				} elsif ($nt1 eq $snp->_nt2 && $nt2 eq $snp->_nt1) {
					${$snps{$pos}->_r_var_info} = 1;
					if ($phased == 1) {
						${$snps{$pos}->_r_phase_info} = -1;
					}
					if ($ps ne "" && exists($PS2blockID{$ps})) {
						${$snps{$pos}->_r_block_id} = $PS2blockID{$ps};
					}
				} else {
					${$snps{$pos}->_r_var_info} = -1;
				}
			} else {
				if ($nt1 eq $snp->_nt1 && $nt2 eq $snp->_nt2) {
					${$snps{$pos}->_r_var_info} = 1;
				} else {
					${$snps{$pos}->_r_var_info} = -1;
				}
			}
		} else {
			# position does not exist in golden set
			inc_count("pred.not_in_ref");
			if ($type eq "het") {
				inc_count("pred.not_in_ref.het");
				if ($phased == 1) {
					inc_count("pred.not_in_ref.het.phased");
				} else {
					inc_count("pred.not_in_ref.het.unphased");
				}
			} elsif ($type eq "hom") {
				inc_count("pred.not_in_ref.hom");
			} elsif ($type eq "ref") {
				inc_count("pred.not_in_ref.ref");
			} else {
				inc_count("pref.not_in_ref.unknown_type");
			}
		}
		if ($type eq "het") {
			inc_count("pred.het");
			if ($phased == 1) {
				inc_count("pred.het.phased");
			} else {
				inc_count("pred.het.unphased");
			}
		} elsif ($type eq "hom") {
			inc_count("pred.hom");
		} elsif ($type eq "ref") {
			inc_count("pred.ref");
		} else {
			inc_count("pref.unknown_type");
		}
		inc_count("pred.total");
	}
	close($fh);

	return $num_PS;
}



##
#	Load golden standard set
##
sub load_golden_set {
	my $file = shift;
	my $fh_r;
	
	if ($file=~m/\.gz$/) {
		$fh_r = new IO::Zlib;
		$fh_r->open("$file", "rb") or die "cannot open file: $!";
	} else {
		open($fh_r, "<$file") or die "cannot open file: $!";
	}
	while (<$fh_r>) {
		if (m/^#/) {
			next;
		}
		chomp;
		my ($chr, $pos, $t1, $ref, $alt, $t2, $t3, $t4, $t5, $gt) = split/\t/, $_;
		if ($chr_check eq "") {
			$chr_check = $chr;
		} elsif ($chr_check ne $chr) {
			print STDERR "multiple chromosome names exist: $_\n";
			exit;
		}
		if (length($ref)>1) {
			print STDERR "Non-SNP variant exists: $_\n";
			exit;
		}
		my @alts = split/,/, $alt;
		for (my $i=0; $i<=$#alts; $i++) {
			if (length($alts[$i])>1) {
				print STDERR "Non-SNP variant exists: $_\n";
				exit;
			}
		}
	
		my $phased = 0;
		my ($nt1, $nt2, $type);
	
		if ($gt =~ m/^(\d+)([\|\/])(\d+)/) {
			my $idx1 = $1;
			my $idx2 = $3;
			if ($2 eq "|") {
				$phased = 1;
			}
			if ($idx1 == $idx2) { # hom
				if ($idx1 == 0) {
					$nt1 = $nt2 = $ref;
					$type = "ref";
				} else {
					$nt1 = $nt2 = $alts[$idx1-1];
					$type = "hom";
				}
			} else { # het
				$type = "het";
				if ($idx1 == 0) {
					$nt1 = $ref;
				} else {
					$nt1 = $alts[$idx1-1];
				}
				if ($idx2 == 0) {
					$nt2 = $ref;
				} else {
					$nt2 = $alts[$idx2-1];
				}
			}
		} else {
			print STDERR "unexpected format: $_\n";
			exit;
		}
	
		my $snp = SNP->new(_pos=>$pos,
							_phased=>$phased,
							_type=>$type,
							_nt1=>$nt1,
							_nt2=>$nt2,
							_nt_ref=>$ref,
							_r_block_id=>0,
							_r_phase_info=>0,
							_r_var_info=>0,
							_r_first_in_block=>0,
							_r_error=>0,
							_r_last_in_block=>0);
		$snps{$pos} = $snp;
	}
	close($fh_r);
}
##


