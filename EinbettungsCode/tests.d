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

void main(){
    auto g = Agraph!2(10,[[1,2],[1,3],[1,7],[2,3],[2,4],[3,7],[4,7],[4,8],[4,10],[5,7],[5,8],[6,8],[6,9],[8,9],[9,10]]);
    auto embeddings = new Embedding(g,true,[],true,true);
    embeddings.check_legal();
    embeddings.logs.writeln;
    embeddings.combinatorial.faces().writeln;
    //embeddings.ship_embedding("wow");
    //return;
    //auto g = Agraph!2([1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16],[[2,5],[2,9],[3,7],[3,16],[4,12],[4,15],[6,8],[7,16],[8,16],[9,12],[10,13],[12,13],[12,15],[13,16]]);
    //auto embeddings = new Embedding(g,true, [],true,true);
    if (!embeddings.check_legal()) {
        writeln("failed - illegal");
        embeddings.logs.writeln;
        embeddings.ship_embedding("failed");
        return;
    } else {
        writeln("legal");
    }
    
    auto graph = Agraph!2([1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17],[[2,5],[2,7],[2,8],[2,12],[2,14],[2,16],[3,7],[3,8],[4,6],[4,8],[4,14],[4,15],[5,14],[6,12],[6,14],[7,12],[7,16],[8,15],[9,11],[9,16],[10,16],[11,16],[12,14],[12,16],[13,15],[14,15],[14,17]]);
    auto embedding = new Embedding(graph,true, [],true,true);
    //auto embedding = new Embedding(graph, [Apair(0,6),Apair(3,6),Apair(3,10),Apair(10,10), Apair(10,0), Apair(-10,0), Apair(-10,10), Apair(-3,10), Apair(-3,6)].reverse);

    //embedding.logs.writeln;
    if (!embedding.check_legal()) {
        writeln("failed - illegal");
        embedding.logs.writeln;
        embedding.ship_embedding("failed");
        return;
    } else {
        writeln("legal");
    }
    graph = Agraph!2([1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20],[[1,7],[1,13],[2,12],[2,13],[3,4],[3,5],[3,10],[3,13],[5,17],[6,19],[8,17],[8,20],[9,20],[10,11],[10,15],[10,19],[11,17],[11,18],[14,20],[17,18]]);
    embedding = new Embedding(graph,true, []);
    //auto embedding = new Embedding(graph, [Apair(0,6),Apair(3,6),Apair(3,10),Apair(10,10), Apair(10,0), Apair(-10,0), Apair(-10,10), Apair(-3,10), Apair(-3,6)].reverse);
    //embedding.logs.writeln;
    if (!embedding.check_legal()) {
        writeln("failed - illegal");
        graph.writeln;
        embedding.logs.writeln;
        embedding.ship_embedding("failed");
        return;
    } else {
        writeln("legal");
    }
    graph = Agraph!2([1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18],[[1,10],[1,13],[2,5],[3,7],[3,12],[4,13],[4,15],[5,14],[8,13],[10,13],[11,13],[11,17],[13,14],[13,15],[14,17],[14,18]]);
    embedding = new Embedding(graph,true, []);
    //auto embedding = new Embedding(graph, [Apair(0,6),Apair(3,6),Apair(3,10),Apair(10,10), Apair(10,0), Apair(-10,0), Apair(-10,10), Apair(-3,10), Apair(-3,6)].reverse);
    //embedding.logs.writeln;
    if (!embedding.check_legal()) {
        writeln("failed - illegal");
        embedding.logs.writeln;
        embedding.ship_embedding("failed");
        return;
    } else {
        writeln("legal");
    }
}