% 03_DTM_unanchored_seed20_mm2.pat: 03_DTM_unanchored_seed20 with up to 2 substitutions in the right TIR (~p2[2,0,0]); TSD exact
%@element Mutator_TIR_transposon Name=DTM Classification=TIR/DTM
p1=9...11   %@ target_site_duplication
p2=20...20  %@ five_prime_terminal_inverted_repeat
50...30000
~p2[2,0,0]  %@ three_prime_terminal_inverted_repeat
p1          %@ target_site_duplication
