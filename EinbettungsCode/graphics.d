mixin(import("aloft.all"));
import std.stdio :writeln, File;
import embedding;
import std.conv;
import std.algorithm.mutation: remove, reverse;
import std.typecons: Tuple,tuple;
import std.math;



void main() {
    auto border = [Apair(0,6),Apair(3,6),Apair(3,10),Apair(10,10), Apair(10,0), Apair(-10,0), Apair(-10,10), Apair(-3,10), Apair(-3,6)].reverse;
    auto path = Agraph!2(8);
    path.insert(1,[2,7]);
    path.insert(2,[3,8]);
    path.insert(3,4);
    path.insert(4,5);
    path.insert(5,6);
    path.insert(6,7);
    path.insert(8,5);
    auto embed = new Embedding(path,true,border);
    embed.ship_embedding("../Arbeit/img/Pfad_eingebettet");
    auto circle = Agraph!2(7);
    circle.insert(1,[2,7]);
    circle.insert(2,3);
    circle.insert(3,4);
    circle.insert(4,5);
    circle.insert(5,6);
    circle.insert(6,7);
    embed = new Embedding(circle,true,border);
    embed.ship_embedding("../Arbeit/img/Kreis_eingebettet");
    auto circle_circle = Agraph!2(9);
    circle_circle.insert(1,[2,7]);
    circle_circle.insert(2,3);
    circle_circle.insert(3,[4,8,9]);
    circle_circle.insert(4,5);
    circle_circle.insert(5,6);
    circle_circle.insert(6,7);
    circle_circle.insert(8,9);
    embed = new Embedding(circle_circle,true,border);
    if (!embed.planar || !embed.check_legal()) {
        embed.logs.writeln;
        return;
    }
    embed.ship_embedding("../Arbeit/img/Kreis_verbunden_eingebettet");
    auto line_demo = Agraph!2(6);
    line_demo.insert(1,[2,3,4]);
    line_demo.insert(2,3);
    line_demo.insert(4,5);
    line_demo.insert(5,6);
    embed = new Embedding(line_demo,true, [],false,false,12.0,true,false);
    embed.ship_embedding("../Arbeit/img/Line_demo");
    auto comp_relview = Agraph!2(6);
    comp_relview.insert(1,[2,5,6]);
    comp_relview.insert(2,[3,4]);
    comp_relview.insert(3,[4,6]);
    comp_relview.insert(4,5);
    comp_relview.insert(5,6);
    embed = new Embedding(comp_relview,true);
    embed.ship_embedding("../Arbeit/img/comp_relview");
    auto comp_konvex = Agraph!2(6);
    comp_konvex.insert(1,[2,3,4,5]);
    comp_konvex.insert(2,[3,6]);
    comp_konvex.insert(3,6);
    comp_konvex.insert(4,[5,6]);
    comp_konvex.insert(5,6);
    embed = new Embedding(comp_konvex,true,[]);
    embed.ship_embedding("../Arbeit/img/comp_konvex");
    auto comp_vinets = Agraph!2(32);
    comp_vinets.insert(1,[2,10,11,25,26]);
    comp_vinets.insert(2,[3,15]);
    comp_vinets.insert(3,[4,16,17]);
    comp_vinets.insert(4,[5,18]);
    comp_vinets.insert(5,[6,19]);
    comp_vinets.insert(6,[7,22]);
    comp_vinets.insert(7,[8,23]);
    comp_vinets.insert(8,[9,23,24]);
    comp_vinets.insert(9,[10,24]);
    comp_vinets.insert(10,[25]);
    comp_vinets.insert(11,[12,23]);
    comp_vinets.insert(12,[13,23]);
    comp_vinets.insert(13,[14,27]);
    comp_vinets.insert(14,[15,28]);
    comp_vinets.insert(15,[16,30]);
    comp_vinets.insert(16,[30,31]);
    comp_vinets.insert(17,[18,31]);
    comp_vinets.insert(18,19);
    comp_vinets.insert(19,20);
    comp_vinets.insert(20,[21,30,31]);
    comp_vinets.insert(21,[22,32]);
    comp_vinets.insert(22,[23,27,29]);
    comp_vinets.insert(23,26);
    comp_vinets.insert(24,25);
    comp_vinets.insert(25,26);
    comp_vinets.insert(27,29);
    comp_vinets.insert(28,[29,30]);
    comp_vinets.insert(29,32);
    comp_vinets.insert(30,32);
    embed = new Embedding(comp_vinets,true,[]);
    if(!embed.check_legal()) {
        embed.logs.writeln;
        return;
    }
    embed.ship_embedding("../Arbeit/img/comp_vinets");

    auto fast_winkelhalbierende = [Apair(0,0), Apair(1,0), Apair(5,1), Apair(5,5), Apair(1,5), Apair(1,1),Apair(0,1)];
    auto Winkelhalbierende = Als(Apair(1,0),Apair(0.8,1));
    draw_output("../Arbeit/img/SchnellWinkel", [tuple(fast_winkelhalbierende,Ablue)],[tuple([Winkelhalbierende],Agreen)]);

    auto close_winkel = [Apair(0,0), Apair(1,0), Apair(5,1), Apair(1,1), Apair(5,5), Apair(0,5)];
    Winkelhalbierende = Als(Apair(1,0),Apair(0,5));
    draw_output("../Arbeit/img/NahWinkel", [tuple(close_winkel,Ablue)],[tuple([Winkelhalbierende],Agreen)]);

    auto spike = [Apair(0,0), Apair(5,0), Apair(5,4), Apair(2,2), Apair(5,5), Apair(0,5)];
    auto skeleton = generate_straight_skeleton(spike, [], [],1.0,false,false, true);
    Apath[] skeleton_graph = [];
    foreach(u,v ; skeleton) {
        skeleton_graph ~= Als(skeleton.val(u),skeleton.val(v));
    }
    draw_output("../Arbeit/img/SkeletonSpike",[tuple(spike,Ablue)],[tuple(skeleton_graph,Amagenta)]);

    auto polygon_to_split = [Apair(0,0), Apair(4,0), Apair(4,4),Apair(4.5,5),Apair(4,4),Apair(0,4)];
    auto split = split_polygon(polygon_to_split,polygon_to_split);
    draw_output("../Arbeit/img/SplitPolygon",[tuple(split[0],Ablue),tuple(split[2],Agreen)]);

    auto border_to_split = [Apair(0,0), Apair(1,1.5), Apair(2,0)];
    border = [Apair(-1,-1), Apair(3,-1), Apair(1,2.5)];
    split = split_polygon(border_to_split,border);
    draw_output("../Arbeit/img/SplitBorder",[tuple(split[0],Ablue)]);

    border = [Apair(0,0),Apair(5,0), Apair(5,5),Apair(0,5)];
    auto part = partition(border,[3,1,1]);
    draw_output("../Arbeit/img/Partition",[tuple(part[0],Ablue),tuple(part[1],Agreen),tuple(part[2],Amagenta)]);

    auto s = 6.0;
    auto demo_graph = Agraph!2(10,[[8,9],[8,10],[9,10],[1,2],[1,4],[2,3],[3,4],[4,5],[4,6],[5,6],[1,7],[3,7]]);
    embed = new Embedding(demo_graph,true,[],true,true,s);
    embed.ship_embedding("../Arbeit/img/demo_graph");

    border = [];
    for(int i = 0; i < demo_graph.order; i++) {
        border ~= Apair(cos(i * 2 * PI /demo_graph.order), sin(i * 2 * PI / demo_graph.order)) * s;
    }

    auto demo_partition = partition(border,[16,6]);
    draw_output("../Arbeit/img/Demopartition",[tuple(demo_partition[0],Ablue),tuple(demo_partition[1],Agreen)]);

    auto offset_1 = find_circle(demo_partition[0],5.0);
    auto offset_2 = find_circle(demo_partition[1],5.0);
    draw_output("../Arbeit/img/demo_offset",[tuple(offset_1,Ablue),tuple(offset_2,Agreen),tuple(demo_partition[0],Ared),tuple(demo_partition[1],Ared)]);

    auto circle_embedded = Agraph!2(4,[[1,2],[1,4],[2,3],[3,4]]);
    embed = new Embedding(circle_embedded,true,demo_partition[0]);
    embed.ship_embedding("../Arbeit/img/demo_circle");

    auto line_embedded = Agraph!2([1,2,3,4,7],[[1,2],[1,4],[2,3],[3,4],[1,7],[3,7]]);
    embed = new Embedding(line_embedded,true,demo_partition[0]);
    embed.ship_embedding("../Arbeit/img/demo_line");

    auto demo_antenne = [Apair(0,0),Apair(5,0),Apair(5,5),Apair(0,5),Apair(0,0),Apair(1,1),Apair(2,2),Apair(3,3),Apair(2,2),Apair(1,1)];
    draw_output("../Arbeit/img/Antenne_Beispiel",[tuple(demo_antenne,Ablue)]);
    ulong[] indices = [4,5];
    Apair[] removed = simplify_remove_illegal(demo_antenne,indices);
    ulong[] relevant_indices = new ulong[](indices.length);
    relevant_indices[] = removed.length;
    for(int i = 0; i < removed.length; i++) {
        for(int j = 0; j < indices.length; j++) {
            if (demo_antenne[indices[j]].eq(removed[i])) {
                relevant_indices[j] = i;
            }
        }
        bool set = true;
        for(int j = 0; j < relevant_indices.length; j++) {
            if (relevant_indices[j] == removed.length) {
                set = false;
            }
        }
        if (set) {
            break;
        }
    }
    auto simp_p1 = simplify_line(removed,[],relevant_indices);
    draw_output("../Arbeit/img/VereinfachungÜberschneidung",[tuple(demo_antenne,Agreen),tuple(simp_p1,Ablue)]);

    auto skeleleton = generate_straight_skeleton(simp_p1, [], [],1.0,false,false, true);
    Apath[] skeleleton_graph = [];
    foreach(u,v ; skeleleton) {
        skeleleton_graph ~= Als(skeleleton.val(u),skeleleton.val(v));
    }
    draw_output("../Arbeit/img/SkeletonSimped",[tuple(demo_antenne,Agreen),tuple(simp_p1,Ablue)],[tuple(skeleleton_graph,Amagenta)]);

    auto simped = simplify(demo_antenne,[],[4,5]);
    auto simplified_skeleton = generate_straight_skeleton(simped, [], [3,5],1.0,false,false, true);
    Apath[] simplified_skeleton_graph = [];
    foreach(u,v ; simplified_skeleton) {
        simplified_skeleton_graph ~= Als(simplified_skeleton.val(u),simplified_skeleton.val(v));
    }
    //auto weighted_simplified_skeleton = generate_straight_skeleton(simped, [], [3,5],25.0,false,false, true);
    //Apath[] weighted_simplified_skeleton_graph = [];
    //foreach(u,v ; weighted_simplified_skeleton) {
    //    weighted_simplified_skeleton_graph ~= Als(weighted_simplified_skeleton.val(u),weighted_simplified_skeleton.val(v));
    //}

    draw_output("../Arbeit/img/Vereinfachung",[tuple(demo_antenne,Agreen),tuple(simped,Ablue)],[tuple(simplified_skeleton_graph,Amagenta)]);

    auto polyogon = [Apair(0,0), Apair(5,0), Apair(5,5), Apair(0,5)];
    auto hole = [Apair(1,2), Apair(1,3), Apair(4,3), Apair(4,2)];
    auto skeleton_with_hole = generate_straight_skeleton(hole,polyogon,[],1.0,false,false,false);
    Apath[] skeleton_with_hole_graph = [];
    foreach(u,v ; skeleton_with_hole) {
        skeleton_with_hole_graph ~= Als(skeleton_with_hole.val(u),skeleton_with_hole.val(v));
    }
    draw_output("../Arbeit/img/SkeletonHole",[tuple(polyogon,Ablue),tuple(hole,Ablue)],[tuple(skeleton_with_hole_graph,Amagenta)]);

}