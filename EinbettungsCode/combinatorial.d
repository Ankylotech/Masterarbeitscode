mixin(import("aloft.all"));
import std.stdio :writeln, File;
import embedding;
import std.algorithm.sorting: sort;
import std.algorithm.mutation: remove, reverse;
import std.algorithm.comparison : min,max;
import std.algorithm.iteration: uniq;
import std.array;
import std.datetime.stopwatch :  StopWatch, AutoStart;
import std.random : Random, uniform;
import std.conv;
import std.csv;
import std.file: write;
import core.thread.osthread;
import core.time;

void main() {
    int node_tests = 100;
    int total = 1000;
    int counter = 0;
    real density = 0.5;
    auto result = File("../Arbeit/statistics/results-combinatorial.csv","w");
    result.write("nodes,total,planar,unplanar\n");
    writeln("testing ", total, " random Graphs each");
    int[] fails = [];
    for(int i = 2;i <= node_tests; i++) {
        ulong unplanar_graphs = 0;
        ulong planar_time = 0;
        ulong unplanar_time = 0;
        ulong planar_graphs = 0;
        while (counter < total) {
            density = 0.5 * (3 * i - 6.0) / ((i * (i-1)) / 2.0);
            auto graph = Agraph!2(i);
            for(int j = 1; j < i; j++) {
                for(int k = j+1; k < i+1; k ++) {
                    if( uniform01(rnd) < density) {
                        graph.insert(j,k);
                    }
                }
            }
            auto sw = StopWatch(AutoStart.no);
            sw.start();
            auto e = new Embedding(graph,false);
            sw.stop();
            auto time = sw.peek.total!"msecs";
            if (e.planar){
                planar_time += time;
                planar_graphs ++;
            } else {
                unplanar_time += time;
                unplanar_graphs++;
            }
            counter++;
        }
        while(planar_graphs < total) {
            density = (3 * i - 6.0) / ((i * (i-1)) / 2.0);
            auto graph = Agraph!2(i);
            for(int j = 1; j < i; j++) {
                for(int k = j+1; k < i+1; k ++) {
                    if(uniform01(rnd) < density) {
                        auto edge = graph.insert(j,k);
                        auto e = new Embedding(graph,false);
                        if (!e.planar) {
                            graph.removeedge(edge);
                        }
                    }
                }
            }
            auto sw = StopWatch(AutoStart.no);
            sw.start();
            auto e = new Embedding(graph,false);
            sw.stop();
            auto time = sw.peek.total!"msecs";
            planar_time += time;
            planar_graphs ++;
        }
        writeln("avg time for Graphs with ",i ," nodes: ", (planar_time + unplanar_time) / real(unplanar_graphs + planar_graphs), "ms");
        if (fails.length > 0) {
            writeln("we failed tests ",fails);
            break;
        } else {
            result.write(to!string(i) ~","~ to!string((planar_time + unplanar_time) / real(unplanar_graphs + planar_graphs)) ~ "," ~ to!string((planar_time) / real(planar_graphs)) ~ "," ~ to!string((unplanar_time) / real(max(unplanar_graphs,1))) ~ "\n");
        }
        counter = 0;
    }
    result.close();
}