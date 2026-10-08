% 09_DTC_CACTG_seed10.pat: labelled copy for --format (scan_for_matches 0.3.0)
%@element CACTA_TIR_transposon Name=DTC Classification=TIR/DTC
p1=3...3               %@ target_site_duplication
CACTG[0,0,0] p2=5...5  %@ five_prime_terminal_inverted_repeat
50...30000
~p2 CAGTG[0,0,0]       %@ three_prime_terminal_inverted_repeat
p1                     %@ target_site_duplication
