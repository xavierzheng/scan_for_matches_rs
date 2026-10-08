% 03_DTM_G_seed20.pat: labelled copy for --format (scan_for_matches 0.3.0)
%@element Mutator_TIR_transposon Name=DTM Classification=TIR/DTM
p1=9...11            %@ target_site_duplication
G[0,0,0] p2=19...19  %@ five_prime_terminal_inverted_repeat
50...30000
~p2 C[0,0,0]         %@ three_prime_terminal_inverted_repeat
p1                   %@ target_site_duplication
