% 09_DTC_CCC_seed12.pat: labelled copy for --format (scan_for_matches 0.3.0)
%@element CACTA_TIR_transposon Name=DTC Classification=TIR/DTC
p1=2...2             %@ target_site_duplication
CCC[0,0,0] p2=9...9  %@ five_prime_terminal_inverted_repeat
50...30000
~p2 GGG[0,0,0]       %@ three_prime_terminal_inverted_repeat
p1                   %@ target_site_duplication
