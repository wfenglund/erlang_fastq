-module(flowerquest).
-export([split_list/3, hive/1, queen/3, worker/2, explore_flowers/2, start/0]).

%% To do:
%% * Automatically spawn workers.
%% * Make the task more computationally expensive and
%% compare run time with a regular program (python and
%% erlang).

% Take a list and return it split it into chunks of a
% chosen size:
split_list([], _, Out_list) ->
	lists:reverse(Out_list);

split_list(In_list, Size, Out_list) when Size > length(In_list) ->
	lists:reverse([In_list|Out_list]);

split_list(In_list, Size, Out_list) ->
	{H, T} = lists:split(Size, In_list),
	split_list(T, Size, [H|Out_list]).

% Loop and distribute tasks to workers and stop workers
% that requests tasks when all tasks are distributed:
hive(Tasks) ->
	receive
		{giveMeWork, Bee_PID} ->
			case Tasks /= [] of
				true ->
					Bee_PID ! {workTask, lists:nth(1, Tasks)},
					hive(lists:nthtail(1, Tasks));
				false ->
					Bee_PID ! {workStop},
					hive(Tasks)
			end;
		_ ->
			hive(Tasks)
	end.

% Loop until she has recieved all data at which point she
% returns the full results and sends a finishing message
% to start()-function:
queen(Data_package, Start_PID, Tasks) ->
	receive
		{flowerInfo, Results} ->
			case (Tasks - 1) > 0 of
				true ->
					queen(Results ++ Data_package, Start_PID, Tasks - 1);
				false ->
					io:fwrite("~p~n", [lists:sort(Results ++ Data_package)]),
					io:fwrite("~p~n", [length(Results ++ Data_package)]),
					Start_PID ! {allTasksDone}
			end;
		_ ->
			queen(Data_package, Start_PID, Tasks)
	end.

% The task that the workers perform:
explore_flowers([], In_results) ->
	In_results;

explore_flowers([H|T], In_results) ->
	Out_results = [H] ++ In_results,
	explore_flowers(T, Out_results).

% loop and request a task from hive and when done send the
% results to the queen:
worker(Hive_PID, Queen_PID) ->
	Hive_PID ! {giveMeWork, self()},
	receive
		{workTask, Task} ->
			Queen_PID ! {flowerInfo, explore_flowers(Task, [])},
			worker(Hive_PID, Queen_PID);
		{workStop} ->
			ok;
		_ ->
			worker(Hive_PID, Queen_PID)
	end.

start() ->
	List1 = lists:seq(1, 250000),
	Tasks = split_list(List1, 1000, []),
	Hive_PID = spawn(flowerquest, hive, [Tasks]),
	Queen_PID = spawn(flowerquest, queen, [[], self(), length(Tasks)]),
	
	spawn(flowerquest, worker, [Hive_PID, Queen_PID]),
	spawn(flowerquest, worker, [Hive_PID, Queen_PID]),
	spawn(flowerquest, worker, [Hive_PID, Queen_PID]),
	spawn(flowerquest, worker, [Hive_PID, Queen_PID]),
	spawn(flowerquest, worker, [Hive_PID, Queen_PID]),
	
	receive
		{allTasksDone} ->
			io:fwrite("~nall tasks are done.")
	end.
