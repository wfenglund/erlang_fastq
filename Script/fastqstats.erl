-module(fastqstats).
-export([fastq_tester/4, supervisor/2, package_manager/3]).
-export([file_reader/2, worker_launcher/2, gzip_to_binary/1]).
-export([unique_seq_finder/2, seq_highlight/2, print_seqs/3]).
-export([make_seq_table/2, start/0]).

% Helper function to fastq_identifier():
fastq_tester([], _, _, Out_list) ->
	Out_list;

fastq_tester([Element|Remainder], Directory, Suffix, Out_list) ->
	RE_result = re:run(Element, Suffix, [{capture, none}]),
	if
		RE_result == match ->
			Ele_path = Directory ++ Element,
			fastq_tester(Remainder, Directory, Suffix, [Ele_path] ++ Out_list);
		RE_result == nomatch ->
			fastq_tester(Remainder, Directory, Suffix, Out_list)
	end.

% Function that collects forward and reverse packages from
% package_manager() and sends them to start() when all both are
% received:
supervisor(Package, Start_PID) ->
	receive
		{Direction, Results} ->
			case (length(Package) + 1) < 2 of
				true ->
					supervisor([{Direction, Results}] ++ Package, Start_PID);
				false ->
					Start_PID ! {allSeqsDone, [{Direction, Results}] ++ Package}
			end;
		_ ->
			supervisor(Package, Start_PID)
	end.

% Function that receives individual parsed forward and
% reverse files and puts them together and sends the
% collected packages to supervisor():
package_manager({Pack_F, Pack_R}, Start_PID, {Tasks_F, Tasks_R}) ->
	receive
		{forward, Results} ->
			case (Tasks_F - 1) > 0 of
				true ->
					package_manager({Results ++ Pack_F, Pack_R},
						     Start_PID,
						     {Tasks_F - 1, Tasks_R});
				false ->
					Start_PID ! {forward, Results ++ Pack_F},
					package_manager({Pack_F, Pack_R},
						       Start_PID,
						       {Tasks_F, Tasks_R})
			end;
		{reverse, Results} ->
			case (Tasks_R - 1) > 0 of
				true ->
					package_manager({Pack_F, Results ++ Pack_R},
						     Start_PID,
						     {Tasks_F, Tasks_R - 1});
				false ->
					Start_PID ! {reverse, Results ++ Pack_R},
					package_manager({Pack_F, Pack_R},
						       Start_PID,
						       {Tasks_F, Tasks_R})
			end;
		_ ->
			package_manager({Pack_F, Pack_R}, Start_PID, {Tasks_F, Tasks_R})
	end.

% Function that takes a file and parses it and sends the
% resulting list of sequences to a given PID:
file_reader({Direction, File}, Pack_PID) ->
	Bin_list = gzip_to_binary(File),
	Short = lists:sublist(Bin_list, 400000),
	Seq_test = fun({X, _}) -> (X == 2 orelse X rem 4 == 2) end,
	Filtered = lists:filter(Seq_test, lists:zip(lists:seq(1, length(Short)), Short)),
	{_, Unzipped} = lists:unzip(Filtered),
	Seq_short = fun(H, T) -> [string:sub_string(binary_to_list(H), 1, 20)|T] end,
	Shortened = lists:foldl(Seq_short, [], Unzipped),
	Pack_PID ! {Direction, Shortened}.

% Function that takes filenames and starts a process for
% each:
worker_launcher({_, []}, _) ->
	ok;

worker_launcher({Direction, [H|T]}, Pack_PID) ->
	spawn(fastqstats, file_reader, [{Direction, H}, Pack_PID]),
	worker_launcher({Direction, T}, Pack_PID).

% Function that reads a gzipped file, extracts it and splits
% the lines into a list:
gzip_to_binary(File_name) ->
	{ok, File} = file:read_file(File_name),
	Inflated_data = zlib:gunzip(File),
	binary:split(Inflated_data, [<<"\n">>], [global]).

% Function that takes a list of sequences and creates a map
% with unique sequences and counts of each:
unique_seq_finder([], In_map) ->
	In_map;

unique_seq_finder([Seq|Remainder], In_map) ->
	case maps:is_key(Seq, In_map) of
		true ->
			Out_map = maps:merge(In_map, #{Seq=>(maps:get(Seq, In_map)+1)});
		false ->
			Out_map = maps:merge(In_map, #{Seq=>1})
	end,
	unique_seq_finder(Remainder, Out_map).

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

% Main function:
start() ->
	% Set specifications:
	File_folder = "../Raw_data_2/",
	Primer_seq = "AAACTCGTGCCAGCCACC",

	% Identify file paths:
	{ok, File_list} = file:list_dir(File_folder),
	Paths_F = fastq_tester(File_list, File_folder, "1.fq.gz", []),
	Paths_R = fastq_tester(File_list, File_folder, "2.fq.gz", []),
	io:fwrite("Identified ~p fastq-files:~n", [length(Paths_F ++ Paths_R)]),
	lists:foreach(fun(X) -> io:fwrite("~p~n", [X]) end, Paths_F ++ Paths_R),	
	io:fwrite("~n"),
	
	% Spawn processes for concurrent parsing of data:
	io:fwrite("Spawning tasks.~n"),
	Task_lengths = {length(Paths_F), length(Paths_R)},
	Super_PID = spawn(fastqstats, supervisor, [[], self()]),
	Pack_PID = spawn(fastqstats, package_manager, [{[], []}, Super_PID, Task_lengths]),
	worker_launcher({forward, Paths_F}, Pack_PID),
	worker_launcher({reverse, Paths_R}, Pack_PID),

	receive
		{allSeqsDone, Package} ->
			io:fwrite("-All tasks are done.~n")
	end,
	
	% Create maps:
	io:fwrite("Creating maps.~n"),
	Sorted = lists:sort(Package),
	[{forward, For},{reverse, Rev}] = Sorted,
	For_map = unique_seq_finder(For, #{}),
	Rev_map = unique_seq_finder(Rev, #{}),
	io:fwrite("-All maps are done.~n"),
	
	% Print result table:
	io:fwrite("~nForward:~n"),
	io:fwrite("__________________________________~n"),
	make_seq_table(For_map, Primer_seq),
	io:fwrite("~nReverse:~n"),
	io:fwrite("__________________________________~n"),
	make_seq_table(Rev_map, Primer_seq).
