% 03_DTM_G_seed20_indel.pat: 03_DTM_G_seed20 with up to 2 substitutions, 1 deletion and 1 insertion in the right TIR (~p2[2,1,1]); TSD exact
%@element Mutator_TIR_transposon Name=DTM Classification=TIR/DTM
p1=9...11            %@ target_site_duplication
G[0,0,0] p2=19...19  %@ five_prime_terminal_inverted_repeat
50...30000
~p2[2,1,1] C[0,0,0]  %@ three_prime_terminal_inverted_repeat
p1                   %@ target_site_duplication
