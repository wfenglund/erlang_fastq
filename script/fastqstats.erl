-module(fastqstats).
-export([fastq_tester/4, fastq_identifier/2, list_printer/1, gzip_to_binary/1, sequence_extractor/3, unique_seq_finder/2, file_looper/2, seq_map_generator/1, seq_highlight/2, print_seqs/3, make_seq_table/2, start/0]).

% Helper function to fastq_identifier():
fastq_tester([], _, _, Out_list) ->
	Out_list;

fastq_tester([Element|Remainder], Directory, Suffix, Out_list) ->
	RE_result = re:run(Element, Suffix, [{capture, none}]),
	if
		RE_result == match ->
			Ele_path = string:concat(Directory, Element),
			fastq_tester(Remainder, Directory, Suffix, lists:append(Out_list, [Ele_path]));
		RE_result == nomatch ->
			fastq_tester(Remainder, Directory, Suffix, Out_list)
	end.

% Function that returns which files are fastq:
fastq_identifier(Directory, Suffix) ->
	File_tuple = file:list_dir(Directory),
	File_list = element(2, File_tuple),
	fastq_tester(File_list, Directory, Suffix, []).

% Function that prints every element of a list:
list_printer([]) ->
	ok;

list_printer([Element|Remainder]) ->
	io:fwrite("~s~n", [Element]),
	list_printer(Remainder).

% Function that reads a gzipped file, extracts it and splits
% the lines into a list:
gzip_to_binary(File_name) ->
	{ok, File} = file:read_file(File_name),
	Inflated_data = zlib:gunzip(File),
	binary:split(Inflated_data, [<<"\n">>], [global]).

% Function that extract every sequence from a fastq-list:
sequence_extractor([], _, In_list) ->
	In_list;

sequence_extractor([Element|Remainder], Counter, In_list) ->
	case Counter == 1 orelse ((Counter - 1) rem 4 == 0) of
		true ->
			Out_list = [Element|In_list];
		false ->
			Out_list = In_list
	end,
	case Counter =< 400000 of
		true ->
			sequence_extractor(Remainder, Counter + 1, Out_list);
		false ->
			sequence_extractor([], Counter, Out_list)
	end.

% Function that takes a list of sequences and creates a map
% with unique sequences and counts of each:
unique_seq_finder([], Inmap) ->
	Inmap;

unique_seq_finder([Element|Remainder], In_map) ->
	Short_seq = string:sub_string(binary_to_list(Element), 1, 20),
	case maps:is_key(Short_seq, In_map) of
		true ->
			Out_map = maps:merge(In_map, #{Short_seq=>(maps:get(Short_seq, In_map)+1)});
		false ->
			Out_map = maps:merge(In_map, #{Short_seq=>1})
	end,
	unique_seq_finder(Remainder, Out_map).

% Helper function for seq_map_generator():
file_looper([], In_list) ->
	In_list;

file_looper([File|Remainder], In_list) ->
	Bin_list = gzip_to_binary(File),
	Out_list = sequence_extractor(Bin_list, 0, []),
	file_looper(Remainder, lists:merge(In_list, Out_list)).

% Function that loads all fastq-files and extracts the
% individual sequences:
seq_map_generator(File_list) ->
	Seq_list = file_looper(File_list, []),
	unique_seq_finder(Seq_list, #{}).

% Function that takes a sequence and a primer and
% returns the sequence with the overlap in blue:
seq_highlight(Seq, []) ->
	Seq;

seq_highlight(Seq, Primer) ->
	case re:run(Seq, Primer, [{capture, none}]) == match of
		true ->
			re:replace(Seq, Primer, "\e[0;34m" ++ Primer ++ "\e[0;37m");
		false ->
			seq_highlight(Seq, lists:droplast(Primer))
	end. 

% Helper function to make_seq_table():
print_seqs([], _, _) ->
	ok;

print_seqs([{Seq, Num}|Remainder], Primer, Counter) ->
	case Counter < 10 of
		true ->
			io:fwrite("~s\t~10B~n", [seq_highlight(Seq, Primer), Num]),
			print_seqs(Remainder, Primer, Counter + 1);
		false ->
			ok
	end.

% Function that takes a map with sequences and counts and
% generates a table:
make_seq_table(Seq_map, Primer) ->
	io:fwrite("Sequence:\t\tFrequency:~n"),
	Seq_list = maps:to_list(Seq_map),
	Seq_list_sort = lists:keysort(2, Seq_list),
	Seq_list_sort_rev = lists:reverse(Seq_list_sort),
	print_seqs(Seq_list_sort_rev, Primer, 0).

start() ->
	File_folder = "../data/",
	Primer_seq = "AAACTCGTGCCAGCCACC",
	Fastq_paths_f = fastq_identifier(File_folder, "1.fq.gz"),
	Fastq_paths_r = fastq_identifier(File_folder, "2.fq.gz"),
	io:fwrite("Identified fastq-files:~n"),
	list_printer(lists:merge(Fastq_paths_f, Fastq_paths_r)),
	io:fwrite("~n"),
	Seq_map_f = seq_map_generator(Fastq_paths_f),
	Seq_map_r = seq_map_generator(Fastq_paths_r),
	io:fwrite("Forward:~n"),
	io:fwrite("__________________________________~n"),
	make_seq_table(Seq_map_f, Primer_seq),
	io:fwrite("~nReverse:~n"),
	io:fwrite("__________________________________~n"),
	make_seq_table(Seq_map_r, Primer_seq).
