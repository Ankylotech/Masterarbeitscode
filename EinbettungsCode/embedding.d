mixin(import("aloft.all"));
import std.stdio : writeln, write;
import std.algorithm: canFind;
import std.math.trigonometry: atan2, sin, cos, asin, tan;
import std.math.constants: PI;
import std.math: isNaN;
import std.math.exponential: pow;
import std.math.algebraic: abs;
import std.typecons: Tuple,tuple;
import std.algorithm.sorting: sort;
import std.algorithm.mutation: remove, reverse, SwapStrategy;
import std.algorithm.comparison : min,max;
import std.algorithm.iteration: uniq, map;
import std.algorithm.searching: countUntil;
import std.conv;
import std.math.rounding: floor,ceil, round;
import std.process: execute;
import std.array: split;
import std.typecons: Nullable,nullable;
import std.array;
import core.math : sqrt;
import core.exception : ArrayIndexError, AssertError;
import linesegmentintersection : lineSegmentLineSegmentIntersection, segmentsIntersect, orientation;


real size = 12.0;
class CombinatorialEmbedding {
    int[] followers;
    int dimension;

    this(int maxVertex) {
        this.followers = new int[](maxVertex*maxVertex);
        this.dimension = maxVertex;
    }

    int getFollower(int from, int to) {
        assert(from * this.dimension + to < this.followers.length);
        return this.followers[from * this.dimension + to];
    }

    void setFollower(int from, int to, int follower){
        assert(from * this.dimension + to < this.followers.length);
        this.followers[from * this.dimension + to] = follower;
    }

    void print() {
        for(int i = 1; i < this.dimension; i++){
            writeln(i, ": ", this.followers[(i*this.dimension + 1)..(i+1)*this.dimension]);
        }
    }

    int[] face(int i, int j) {
        int[] face = [i];
        int tail = i,head = j;
        while(head != i || this.getFollower(tail,head) != j) {
            if (this.getFollower(tail,head) >= this.dimension || face.length >= this.followers.length) {
                this.print();
                writeln(face);
                assert(0);
            }
            face ~= head;
            int next = this.getFollower(tail,head);
            tail = head;
            head = next;
        }
        return face;
    }

    int[][] faces() {
        //this.print();
        ulong l = this.dimension;
        bool[][] visited = new bool[][](l,l);
        int[][] result = [];
        for(int i = 1; i < l; i++) {
            for(int j = 1; j < l; j++) {
                if (this.getFollower(i,j) >= this.dimension) {
                    assert(0);
                } if (this.getFollower(i,j) != 0 && !visited[i][j]) {
                    auto c = face(i,j);
                    for(int k = 0; k < c.length; k++) {
                        visited[c[k]][c[(k+1) % c.length]] = true;
                    }
                    result ~= c;
                }
            }
        }
        return result;
    }
}

class Embedding {
    CombinatorialEmbedding combinatorial;
    Agnodearray!(Apair) nodepositions;
    Aedgearray!(Apath) edgepositions;
    const Agraph!() G; // The Graph to embed
    Agraph!() H; // The part that is already embedded
    bool planar;
    bool legal;
    bool geometric;
    bool skeleton;
    bool simplified;
    Appender!(string) logs;
    Apair[] border;

    this(const Agraph!() G, bool geometric, Apair[] border = [], bool simplified = false, bool skeleton = false, real s = size, bool weighted = true, bool sorted = true) {
        G.writeln;
        this.combinatorial = new CombinatorialEmbedding(G.vsup);
        this.nodepositions = G.nodearray!Apair;
        this.edgepositions = G.edgearray!Apath;
        this.G = G;
        this.H = Agraph!2();
        this.planar = true;
        this.legal = true;
        this.geometric = geometric;
        this.simplified = simplified;
        this.skeleton = skeleton;
        if (border.length == 0 && this.geometric){
            auto circle_size = max(this.G.order,3);
            for(int i = 0; i < circle_size; i++) {
                border ~= Apair(cos(i * 2 * PI /circle_size), sin(i * 2 * PI / circle_size)) * s;
            }
        }
        this.border = border;
        Agraph!()[] components = [];
        auto visited = new bool[](G.vsup);
        visited[] = false;
        int[] dfs(int current, int last) {
            visited[current] = true;
            int[] result = [current];
            this.G.iter(current, (n) {
                if (n != last && !visited[n]) {
                    result ~= dfs(n,current);
                }
            });
            return result;
        }
        int[] sizes = [];
        foreach(v ; this.G) {
            if (!visited[v]) {
                auto g = this.G[dfs(v,0)];
                sizes ~= g.order + g.nedges;
                components ~= g;
            }
        }
        Apair[][] borders;
        if (this.geometric){
            sort!("a.nedges + a.order > b.nedges + b.order")(components);
            sizes = sizes.sort.array.reverse;
            borders = partition(border,sizes);
        } else {
            borders = new Apair[][](sizes.length);
            borders[] = [];
        }
        foreach(j,g ; components) {
            if (g.order == 1) {
                if (this.geometric){
                    Apair center = Apair(0,0);
                    foreach(p ; borders[j]) {
                        center += p;
                    }
                    center /= borders[j].length;
                    this.setNodePosition(g.choice(), center);
                }
            } else {
                auto bcc_components = g.Abccpartition();
                if (sorted && this.geometric){
                    sort!("a.length > b.length")(bcc_components);
                } else {
                    bcc_components = bcc_components.reverse;
                }
                Agraph!()[] component_graphs = [];
                foreach(bcc ; bcc_components) {
                    int[] bcc_nodes = [];
                    foreach(e ; bcc) {
                        bcc_nodes ~= g.head(e);
                        bcc_nodes ~= g.tail(e);
                    }
                    bcc_nodes = bcc_nodes.sort.uniq.array;
                    if (bcc_nodes.length >= 3 && bcc.length > 3 * bcc_nodes.length - 6) {
                        this.planar = false;
                        return;
                    }
                    component_graphs ~= g[bcc_nodes,bcc];
                }
                Tuple!(ulong,ulong,ulong)[][] connections = calculate_connections(component_graphs);
                ulong[] finished = [];
                ulong[] relevant = [0];
                while(relevant.length > 0 && this.planar) {
                    ulong[] other_indices = [];
                    if (finished.canFind(relevant[0])) {
                        relevant = relevant[1..$];
                        continue;
                    }
                    if (bcc_components[relevant[0]].length == 1) {
                        auto e = bcc_components[relevant[0]][0];
                        auto line = [g.head(e),g.tail(e)];
                        ulong[] connect = [];
                        foreach(l ; connections[relevant[0]]) {
                            connect ~= l[0];
                        }
                        while(connect.length > 0) {
                            auto e2 = bcc_components[connect[0]][0];
                            if (bcc_components[connect[0]].length > 1 || (g.head(e2) != line[0] && g.head(e2) != line[$-1]
                                && g.tail(e2) != line[0] && g.tail(e2) != line[$-1]) 
                                || ((g.head(e2) == line[0] || g.tail(e2) == line[0]) && this.H.isnode(line[0]))
                                || ((g.head(e2) == line[$-1] || g.tail(e2) == line[$-1]) && this.H.isnode(line[$-1]))) {
                                other_indices ~= connect[0];
                            } else {
                                if (g.head(e2) == line[0]) {
                                    line = g.tail(e2) ~ line;
                                } else if (g.tail(e2) == line[0]) {
                                    line = g.head(e2) ~ line;
                                } else if (g.head(e2) == line[$-1]) {
                                    line ~= g.tail(e2);
                                } else {
                                    line ~= g.head(e2);
                                }
                                foreach(conn ; connections[connect[0]]) {
                                    connect ~= conn[0];
                                }
                                finished ~= connect[0];
                            }
                            connect = connect[1..$];
                        }
                        auto face = this.embed_line_combinatorial(line,borders[j]);
                        if(this.geometric) { this.embed_line(line, face, borders[j]);}
                    } else {
                        auto connect = connections[relevant[0]];
                        Tuple!(ulong,ulong)[] embed_data = [];
                        foreach(c ; connect) {
                            if (c[0] != relevant[0]) {
                                other_indices ~= c[0];
                                embed_data ~= tuple(c[1],c[2]);
                            }
                        }
                        if (!this.double_connected_embedding(component_graphs[relevant[0]], borders[j],embed_data, weighted)) {
                            this.planar = false;
                        }
                    }
                    finished ~= relevant[0];
                    relevant = relevant[1..$] ~ other_indices;
                    relevant = relevant.sort.uniq.array;
                }
            }
        }
    }

    Apair[] positions(int[] face) {
        Apair[] result = [];
        foreach(i,v; face) {
            auto p = this.getNodePosition(v);
            Apair last;
            if (result.length == 0 || (!p.eq(result[$-1]))){
                result ~= p;
                last = p;
            } else {
                last = result[$-1];
            }
            auto e = this.G.e(face[i],face[(i+1) % face.length]);
            auto line = this.getEdgePath(e);
            auto total = line.arctime(line.arclength);
            if(!line(0).eq(p)) {
                for(int j = to!int(ceil(total)) - 1; j >= 1; j--) {
                    p = line(j);
                    if (!p.eq(last)) {
                        result ~= p;
                    }
                }
            } else {
                for(int j = 1; j < total; j++) {
                    p = line(j);
                    if (!p.eq(last)) {
                        result ~= p;
                    }
                }
            }
        }
        return result;
    }

    int[] embed_line_combinatorial(int[] current, Apair[] border) {
        int articulation = 0;
        if (this.H.isnode(current[0])) {
            articulation = current[0];
        }
        if (this.H.isnode(current[$-1])) {
            articulation = current[$-1];
            current = current.reverse;
        }
        int[] face = [];
        auto border_area = polygon_area(border);
        if (articulation > 0) {
            auto area = -1.0;
            foreach(n1 ; this.H.neighbors(articulation)) {
                auto f = this.combinatorial.face(articulation,n1);
                if (!this.geometric) {
                    face = f;
                    break;
                } 
                auto a = polygon_area(this.positions(f));
                if (clockwise(this.positions(f))) {
                    a = border_area - a;
                }
                if (a > area) {
                    face = f;
                    area = a;
                }
            }
            this.setFollower(current[1],current[0], face[1]);
            this.setFollower(face[$-1],current[0],current[1]);
            assert(this.H.isnode(current[0]));
        } else {
            this.setFollower(current[1],current[0], current[1]);
            this.H |= current[0];
        }
        this.H |= current[1];
        this.H.insert(current[0],current[1]);
        for(int i = 2; i < current.length; i++) {
            this.H |= current[i];
            this.setFollower(current[i-2],current[i-1],current[i]);
            this.setFollower(current[i],current[i-1],current[i-2]);
            this.H.insert(current[i-1],current[i]);
        }
        this.setFollower(current[$-2],current[$-1],current[$-2]);
        return face;

    }

    void embed_line(int[] current,int[] face, Apair[] border) {
        Apair start; 
        Apair line;
        if (face.length > 0) {
            auto polygon = this.positions(face);
            if (this.skeleton && this.simplified){
                writeln(counter);
                writeln(face);
                debug_output([tuple(polygon,Ayellow)]);
            }            
            auto intersect_point = first_intersection(0, polygon, border)[2];
            intersect_point *= 0.5;
            start = polygon[0];
            line = intersect_point;
        } else {
            auto max_len = 0.0;
            for(int i = 0; i < border.length; i++) {
                auto intersect = first_intersection(i,border,border)[2] * 0.5;
                if (intersect.abs() > max_len) {
                    start = border[i] + intersect * 0.5;
                    line = intersect;
                    max_len = intersect.abs();
                }
            }
            this.setNodePosition(current[0],start);
        }
        auto prev = start;
        for(int i = 1; i < current.length; i++) {
            auto pos = start + i * line / (current.length-1);
            if(pos.isnan()) {
                this.logs.writeln;
                current.writeln;
                face.writeln;
                border.writeln;
                start.writeln;
                line.writeln;
            }
            assert(!pos.isnan());
            this.setNodePosition(current[i],pos);
            this.setEdgePath(this.G.e(current[i-1],current[i]), Als(prev,pos));
            prev = pos;
        }
    }

    bool double_connected_embedding(Agraph!() unembedded, Apair[] border, Tuple!(ulong,ulong)[] connection_data, bool weighted) {
        auto H = Agraph!2();
        int articulation = 0;
        foreach(v ; unembedded) {
            if (this.H.isnode(v)) {
                articulation = v;
                break;
            }
        }
        
        Agraph!()[] fragments(Agraph!() base, Agraph!() fragment_graph) {
            Agraph!()[] result = [];
            auto finished = new bool[](fragment_graph.vsup);
            finished[] = false;
            void fragment(int current, int last, Agraph!() cur_fragment) {
                cur_fragment |= current;
                if (last != 0) {
                    cur_fragment.insert(last,current);
                }
                if (base.isnode(current)) {
                    finished[current] = true;
                    return;
                }
                fragment_graph.iter(current, (n) {
                    if (n != last) {
                        if (!cur_fragment.isnode(n)) {
                            fragment(n,current,cur_fragment);
                        } else if (finished[n]) {
                            cur_fragment.insert(current,n);
                        }
                    }
                });
                finished[current] = true;
            }
            foreach(v ; fragment_graph) {
                if (base.isnode(v)) {
                    fragment_graph.iter(v, (n) {
                        if (base.isnode(n) && base.degree(n) > 0 && v < n) {
                            result ~= Agraph!2([v,n],[[v,n]]);
                        }
                    });
                } else if (!finished[v]) {
                    auto r = Agraph!2();
                    fragment(v,0,r);
                    result ~= r;
                }
            }
            return result;
        }

        int[] find_path(Agraph!() fragment, int[] anchors) {
            auto visited = new bool[](fragment.vsup);
            visited[] = false;
            for(int i = 0; i < anchors.length; i++) {
                visited[anchors[i]] = true;
            }
            int[] dfs(int current, int last, int target) {
                visited[current] = true;
                int[] result = [];
                fragment.iter(current, (v) {
                    if (v != last){
                        if (v == target) {
                            result = [current,target];
                        }
                        if (!visited[v]) {
                            auto r = dfs(v,current, target);
                            if (r.length > 0 && r.length + 1 > result.length) {
                                result = current ~ r;
                            }
                        }
                        
                    }
                });
                return result;
            }
            return dfs(anchors[0],0,anchors[1]);
        }

        int[] find_cycle() {
            auto visited = new bool[](unembedded.vsup);
            visited[] = false;
            int[] dfs(int current, int last, int target) {
                visited[current] = true;
                int[] result = [];
                unembedded.iter(current, (v) {
                    if (v != last){
                        if (v == target && result.length == 0) {
                            result = [current];
                        } else if (!visited[v]) {
                            auto r = dfs(v,current, target);
                            if (r.length > 0 && r.length + 1 > result.length) {
                                result = current ~ r;
                            }
                        }
                        
                    }
                });
                return result;
            }
            int s = unembedded.choice();
            if (articulation != 0) {
                s = articulation;
            } 
            return dfs(s,0, s);
        }

        int[] anchors(int[] cycle, Agraph!() fragment) {
            int[] result = [];
            foreach(v ; cycle) {
                if (fragment.isnode(v) && fragment.degree(v) > 0) {
                    result ~= v;
                }
            }
            result = result.sort.uniq.array;
            return result;
        }
        auto count = unembedded.nedges;
        if (count > 3 * unembedded.order - 6) { 
            return false;
        }
        auto cycle = find_cycle();
        assert(cycle.length > 2);
        auto border_area = polygon_area(border);
        this.embed_cycle_combinatorial(H, unembedded, cycle, articulation, border_area);
        if (this.geometric) {this.embed_cycle(cycle, border, articulation);}
        if (H.nedges >= count) {
            return true;
        }
        struct FragmentData {
            Agraph!() fragment;
            int[][] fully_connected_faces;
            int[] anchors;
        }
        auto fragment_array = fragments(H,unembedded);
        FragmentData[] B = [];
        foreach(f ; fragment_array) {
            FragmentData fd;
            fd.fragment = f;
            fd.anchors = anchors(cycle,f);
            fd.fully_connected_faces = [cycle, this.combinatorial.face(cycle[1],cycle[0])];
            B ~= fd;
        }
        while(H.nedges < count) {
            bool one_max = false;
            ulong max_edges = 0;
            ulong fragment_index = 0;
            int[][] faces = null;
            foreach(i,f; B) {
                if (one_max && f.fully_connected_faces.length > 1) {
                    continue;
                }
                if(f.fragment.nedges > max_edges || (f.fully_connected_faces.length == 1 && !one_max)) {
                    faces = f.fully_connected_faces;
                    max_edges = f.fragment.nedges;
                    fragment_index = i;
                    if (f.fully_connected_faces.length == 1) {
                        one_max = true;
                    }
                }
            }
            if (faces.length == 0) {
                return false;
            }
            int[] path = find_path(B[fragment_index].fragment, B[fragment_index].anchors);
            int[] face = faces[0];
            Apair[] polygon;
            bool outer;
            if (this.geometric){
                polygon = this.positions(face);
                outer = false;
                real max = polygon_area(polygon);
                if(clockwise(polygon)) {
                    outer = true;
                    max = border_area - max;
                } 
                for(int i = 1; i < faces.length; i++) {
                    int[] next = faces[i];
                    Apair[] poly = this.positions(next);
                    real a = polygon_area(poly);
                    auto _outer = false;
                    if (clockwise(poly)) {
                        _outer = true;
                        a = border_area - a;
                    }
                    if(a > max) {
                        face = next;
                        polygon = poly;
                        outer = _outer;
                        max = a;
                    }
                }
            }
            this.embed_path_combinatorial(H, unembedded, path, face);
            Agraph!() path_g = Agraph!2();
            foreach(i,v ; path) {
                path_g |= v;
                if (i > 0) {
                    path_g.insert(v,path[i-1]);
                    B[fragment_index].fragment.removeedge(v,path[i-1]);
                }
            }
            foreach(i,v ; face) {
                path_g |= v;
                if (i > 0) {
                    path_g.insert(v, face[i-1]);
                }
            }
            path_g.insert(face[0],face[$-1]);
            auto new_fragments = fragments(path_g,B[fragment_index].fragment);
            ulong first_index = face.length, last_index = face.length;
            for (int i = 0; i < face.length; i++) {
                if (face[i] == path[0]) {
                    first_index = i;
                    if (last_index < face.length) {
                        break;
                    }
                }else if (face[i] == path[$-1]) {
                    last_index = i;
                    if (first_index < face.length) {
                        break;
                    }
                }
            }
            int[][] split_faces = [[],[]];
            int k = 0;
            for(int i = 0; i < face.length; i++) {
                if (i == first_index) {
                    split_faces[k] ~= path.dup;
                    k = (k+1) % 2;
                } else if (i == last_index) {
                    split_faces[k] ~= path.dup.reverse;
                    k = (k+1) % 2;
                } else {
                    split_faces[k] ~= face[i];
                }
            }
            B = B.remove!(SwapStrategy.unstable)(fragment_index);
            for(int l = 0; l < B.length; l++) {
                for(int i = 0; i < B[l].fully_connected_faces.length; i++) {
                    if (eq_shift(B[l].fully_connected_faces[i],face)) {
                        B[l].fully_connected_faces = B[l].fully_connected_faces.remove!(SwapStrategy.unstable)(i);
                        bool f1 = true;
                        bool f2 = true;
                        foreach(v ; B[l].anchors) {
                            if (!split_faces[0].canFind(v)) {
                                f1 = false;
                            }
                            if (!split_faces[1].canFind(v)) {
                                f2 = false;
                            }
                            if (!f1 && !f2) {break;}
                        }
                        if (f1) {
                            B[l].fully_connected_faces ~= split_faces[0];
                        }
                        if (f2) {
                            B[l].fully_connected_faces ~= split_faces[1];
                        }
                        break;
                    }
                }
                if (B[l].fully_connected_faces.length == 0) {
                    return false;
                }
            }
            foreach (f ; new_fragments) {
                FragmentData fd;
                fd.fragment = f;
                fd.anchors = anchors(face ~ path[1..$-1],f);
                fd.fully_connected_faces = [];
                bool f1 = true;
                bool f2 = true;
                foreach(v ; fd.anchors) {
                    if (!split_faces[0].canFind(v)) {
                        f1 = false;
                    }
                    if (!split_faces[1].canFind(v)) {
                        f2 = false;
                    }
                    if (!f1 && !f2) {break;}
                }
                if (f1) {
                    fd.fully_connected_faces ~= split_faces[0];
                }
                if (f2) {
                    fd.fully_connected_faces ~= split_faces[1];
                }
                B ~= fd;
            }
            ulong c0 = 0;
            ulong c1 = 0;
            for(int l = 0; l < B.length; l++) {
                for(int i = 0; i < B[l].fully_connected_faces.length; i++) {
                    if (eq_shift(B[l].fully_connected_faces[i],split_faces[0])) {
                        c0 += B[l].fragment.nedges;
                    }
                    if (eq_shift(B[l].fully_connected_faces[i],split_faces[1])) {
                        c1 += B[l].fragment.nedges;
                    }
                }
            }
            if (this.geometric) {
                foreach(c ; connection_data) {
                    if (split_faces[0].canFind(c[0])) {
                        c0 += c[1];
                    } 
                    if (split_faces[1].canFind(c[0])) {
                        c1 += c[1];
                    } 
                }
                this.embed_path(path, face, first_index, last_index, polygon, outer ? border : [], cast(float)(c0 + path.length) / cast(float)(c1 + path.length), weighted);
            }
        }
        return true;
    }

    int getFollower(int from, int to) {
        return this.combinatorial.getFollower(from,to);
    }

    void setFollower(int from, int to, int follower){
        this.combinatorial.setFollower(from,to,follower);
    }

    void setNodePosition(int node, Apair position) {
        this.nodepositions[node] = position;
    }

    Apair getNodePosition(int node) {
        if (!this.geometric) {
            return Apair(0,0);
        }
        return this.nodepositions[node];
    }

    void setEdgePath(int edge, Apath path) {
        if (edge == 0) {
            assert(false);
        }
        assert(edge != 0);
        this.edgepositions[edge] = path;
    }

    Apath getEdgePath(int edge) {
        if (!this.geometric) {
            return Als(Apair(0,0),Apair(0,0));
        }
        return this.edgepositions[edge];
    }

    Agnodearray!(Apair) getNodes() {
        return this.nodepositions;
    }

    Aedgearray!(Apath) getEdges() {
        return this.edgepositions;
    }

    void embed_cycle_combinatorial(Agraph!() H, Agraph!() unembedded, int[] cycle, int articulation, real border_area) {
        int neighbour = 0;
        long articulation_index = -1;
        if (articulation > 0) {
            articulation_index = cycle.countUntil(articulation);
            if (this.geometric && false) {
                auto max_area = 0.0;
                foreach(n ; this.H.neighbors(articulation)) {
                    auto face = this.combinatorial.face(articulation,n);
                    auto polygon = this.positions(face);
                    auto area = polygon_area(polygon);
                    if (clockwise(polygon)) {
                        area = border_area - area;
                    }
                    if (area > max_area) {
                        max_area = area;
                        neighbour = n;
                    }
                }
            } else {
                neighbour = this.H.neighbors(articulation)[0];
            }
        }
        H |= cycle[0];
        this.H |= cycle[0];
        for(int i = 0; i < cycle.length; i++) {
            this.setFollower(cycle[i], cycle[(i+1)%cycle.length], cycle[(i+2)%cycle.length]);
            this.setFollower(cycle[(i+2)%cycle.length], cycle[(i+1)%cycle.length], cycle[i]);
            H |= cycle[(i+1) % cycle.length];
            this.H |= cycle[(i+1) % cycle.length];
            H.insert(cycle[i], cycle[(i+1)%cycle.length]);
            this.H.insert(cycle[i], cycle[(i+1)%cycle.length]);
            unembedded.removeedge(cycle[i], cycle[(i+1) % cycle.length]);
        }
        if (articulation > 0) {
            int other = this.getFollower(neighbour, cycle[articulation_index]);
            this.setFollower(neighbour,cycle[articulation_index], cycle[(articulation_index - 1 + cycle.length) % cycle.length]);
            this.setFollower(cycle[(articulation_index + 1) % cycle.length],cycle[articulation_index], other);
        }
    }
    
    void embed_cycle(int[] cycle, Apair[] border, int articulation) {
        Apair[] circle;
        ulong articulation_index = 0;
        if (articulation > 0) {
            articulation_index = cycle.countUntil(articulation);
            auto neighbour = this.getFollower(cycle[(articulation_index + 1) % cycle.length], articulation);
            int[] face = this.combinatorial.face(articulation,neighbour);
            ulong last = 0;
            for (ulong i = face.length; i > 0; i--) {
                if (face[i - 1] == articulation) {
                    last = i - 1;
                    break;
                }
            }
            face = face[0..last];
            auto polygon = this.positions(face);
            circle = find_circle_with_articulation(polygon,border, this.simplified, this.skeleton);
        } else {
            circle = find_circle(border,5.0, this.simplified, this.skeleton);
        }
        auto p = Als(circle[0],circle[1]);
        for(int i = 1; i < circle.length; i++) {
            p ~= Als(circle[i],circle[(i + 1) % circle.length]);
        }
        if (this.skeleton && this.simplified) {
            writeln(cycle);
        }
        real[] positions = place_along_path(p,cycle.length + 1);
        for(ulong i = 0; i < cycle.length; i++) {
            auto index = (i + articulation_index) % cycle.length; 
            this.setNodePosition(cycle[index],p(positions[i]));
            auto e = this.G.e(cycle[index], cycle[(index+1) % cycle.length]);
            if(floor(positions[i+1]) - floor(positions[i]) >= 1) {
                int c = (cast(int) floor(positions[i])) + 1;
                auto l = Als(p(positions[i]),p(c));
                while(c+1 < positions[i+1]) {
                    l ~= Als(p(c),p(c+1));
                    c += 1;
                }
                l ~= Als(p(c),p(positions[i+1]));
                this.setEdgePath(e,l);
            } else {
                this.setEdgePath(e,Als(p(positions[i]),p(positions[i+1])));
            }
        }

    }

    void embed_path_combinatorial(Agraph!() H, Agraph!() unembedded, const int[] path, const int[] face) {
        ulong first_index = face.length, last_index = face.length;
        for (int i = 0; i < face.length; i++) {
            if (face[i] == path[0]) {
                first_index = i;
                if (last_index < face.length) {
                    break;
                }
            }else if (face[i] == path[$-1]) {
                last_index = i;
                if (first_index < face.length) {
                    break;
                }
            }
        }
        
        for(int i = 1;i + 1 < path.length; i++) {
            H |= path[i];
            H.insert(path[i-1],path[i]);
            this.H |= path[i];
            this.H.insert(path[i-1],path[i]);
            unembedded.removeedge(path[i-1],path[i]);

            this.setFollower(path[i-1],path[i],path[i+1]);
            this.setFollower(path[i+1],path[i],path[i-1]);
        }
        H.insert(path[$-2],path[$-1]);
        this.H.insert(path[$-2],path[$-1]);
        unembedded.removeedge(path[$-2],path[$-1]);
        this.setFollower(face[(first_index + face.length-1) % face.length], path[0], path[1]);
        this.setFollower(path[1],path[0], face[(first_index + 1) % face.length]);
        this.setFollower(path[$-2],path[$-1], face[(last_index + 1) % face.length]);
        this.setFollower(face[(last_index + face.length - 1) % face.length], path[$-1], path[$-2]);
    }

    void embed_path(const int[] path, const int[] face, ulong first_index, ulong last_index,const Apair[] face_polygon, const Apair[] outer, real factor, bool weighted) {
        assert(clockwise(face_polygon.dup) == (outer.length != 0));
        ulong ff_index = face_polygon.length, lf_index = face_polygon.length;
        for(int i = 0; i < face_polygon.length; i++) {
            if(face_polygon[i] == this.getNodePosition(face[first_index])) {
                if (ff_index >= face_polygon.length || lf_index >= face_polygon.length) {
                    ff_index = i;
                }
            } else if (face_polygon[i] == this.getNodePosition(face[last_index])) {
                if (ff_index >= face_polygon.length || lf_index >= face_polygon.length) {
                    lf_index = i;
                }
            }
        }
        auto p = find_polygon_path(face_polygon, outer, ff_index, lf_index, factor, this.simplified, this.skeleton, weighted);
        assert(!p[0].isnan());
        assert(!p[1].isnan());
        auto pp = Als(p[0],p[1]);
        for(int i = 1;i + 1 < p.length; i++) {
            assert(!p[i+1].isnan());
            pp ~= Als(p[i],p[i+1]);
        }
        real[] times = place_along_path(pp,path.length);
        auto prev = 0.0;
        auto prev_pos = this.getNodePosition(path[0]);
        for(int i = 1;i + 1 < path.length; i++) {
            auto e = this.G.e(path[i-1],path[i]);
            auto next = pp(times[i]);
            this.setNodePosition(path[i],next);
            if(floor(times[i]) - floor(prev) >= 1) {
                int c = (cast(int) floor(prev)) + 1;
                auto l = Als(prev_pos,pp(c));
                while(c+1 < times[i]) {
                    l ~= Als(pp(c),pp(c+1));
                    c += 1;
                }
                l ~= Als(pp(c),pp(times[i]));
                this.setEdgePath(e,l);
            } else {
                this.setEdgePath(e,Als(prev_pos,next));
            }
            prev = times[i];
            prev_pos = next;
        }
        auto e = this.G.e(path[$-2],path[$-1]);
        if(times[$-1] - floor(prev) >= 1) {
            int c = (cast(int) floor(prev)) + 1;
            auto l = Als(prev_pos,pp(c));
            while(c+1 < times[$-1]) {
                l ~= Als(pp(c),pp(c+1));
                c += 1;
            }
            l ~= Als(pp(c),this.getNodePosition(path[$-1]));
            this.setEdgePath(e,l);
        } else {
            this.setEdgePath(e,Als(prev_pos,this.getNodePosition(path[$-1])));
        }
    }

    void ship_embedding(string name) {
        if (!this.geometric) {return;}
        auto finale = Agembedding!()(this.getNodes(),this.getEdges());
        draw(this.G.visual!(1,0)(finale));
        draw(this.border,Ared);
        for(int i = 0; i < this.border.length; i++) {
            draw(Als(border[i],border[(i+1) % border.length]), Ared);
        }
        shiperase(name);
    }

    bool check_legal() {
        if (!this.legal) {
            return false;
        }
        if (!this.geometric || !this.planar) {
            return true;
        }
        auto edges = this.G.edgeset;
        if (edges.length == 0) {
            return true;
        }
        foreach(e1 ; edges) {
            auto h1 = this.G.head(e1);
            auto t1 = this.G.tail(e1);
            auto p1 = this.getEdgePath(e1);
            foreach(e2 ; edges) {
                if (e1 == e2) {
                    continue;
                }
                auto h2 = this.G.head(e2);
                auto t2 = this.G.tail(e2);
                auto p2 = this.getEdgePath(e2);
                auto intersections = intersection_paths(p1,p2);
                foreach(intersect; intersections){
                    if (!(this.getNodePosition(h1).eq(intersect,Aeps.e[3],0.01)) && !(this.getNodePosition(t1).eq(intersect,Aeps.e[3],0.01))) {
                        this.logs ~= "collision between line " ~ to!string(h1) ~ "," ~ to!string(t1) ~ " and " ~ to!string(h2) ~ "," ~ to!string(t2) ~ " and position " ~ to!string(intersect) ~ "\n";
                        //this.logs ~= to!string(p1.captimes(p2)) ~ " " ~ to!string(p1) ~ " " ~ to!string(p2) ~ "\n";
                        this.logs ~= "distance to ends: " ~ to!string((this.getNodePosition(h1) - intersect).abs()) ~ " " ~ to!string((this.getNodePosition(t1) - intersect).abs()) ~ "\n";
                        this.logs ~= to!string(this.getNodePosition(h1)) ~ "."  ~ to!string(this.getNodePosition(t1)) ~ "\n";
                        this.logs ~= to!string(intersections) ~ "\n";
                        this.logs ~= to!string(p1) ~ " " ~ to!string(p2) ~ "\n";
                        return false;
                    }
                }
                
            }
            foreach(v ; this.G) {
                if (h1 == v || t1 == v) {
                    continue;
                }
                auto p = this.getNodePosition(v);
                assert(!p.isnan());
                if (p1.closestpoint(p) == p) {
                    this.logs ~= "The line " ~ to!string(h1) ~ "," ~ to!string(t1) ~ " goes through the point " ~ to!string(v) ~ " at " ~ to!string(p) ~ "\n";
                    return false;
                }
            }
        }
        return true;
    }
}
Tuple!(ulong,ulong,ulong)[][] calculate_connections(Agraph!()[] graphs) {
    Tuple!(ulong,ulong,ulong)[][] result = new Tuple!(ulong,ulong,ulong)[][](graphs.length);
    result[] = [];
    bool[] visited = new bool[](graphs.length);
    ulong[] current = [0];
    while(current.length > 0) {
        auto c = current[$-1];
        auto graph = graphs[c];
        ulong[] added = [];
        foreach(v ; result[c]){
            added ~= v[0];
        }
        bool visit = true;
        for(ulong i = 0; i < graphs.length; i++) {
            if (i == c || added.canFind(i)) {
                continue;
            }
            foreach (v ; graph)
            {
                if (graphs[i].isnode(v)) {
                    if (visited[i]) {
                        ulong g_sum = 0;
                        foreach(r ; result[i]) {
                            g_sum += r[2];
                        }
                        g_sum += graphs[i].order + graphs[i].nedges;
                        result[c] ~= tuple(i,to!ulong(v),g_sum);
                    } else if (!current.canFind(i)){
                        current ~= i;
                        visit = false;
                    }
                    break;
                }
            }
        }
        if (visit) {
            visited[c] = visit;
            current = current[0..$-1];
        }
    }
    return result;
}

bool eq_shift(int[] f1, int[] f2) {
    if (f1.length != f2.length) {
        return false;
    }
    int m1 = f1[0];
    int i1 = 0;
    int m2 = f2[0];
    int[] i2 = [0];
    for(int i = 1; i < f1.length; i++) {
        if (f1[i] < m1) {
            m1 = f1[i];
            i1 = i;
        } 
        if (f2[i] < m2) {
            m2 = f2[i];
            i2 = [i];
        } else if (f2[i] == m2) {
            i2 ~= i;
        }
    }
    if (m1 != m2) {
        return false;
    }
    foreach(j ; i2){
        bool match = true;
        for(int i = 0; i < f1.length; i++) {
            if (f1[i] != f2[(i + f2.length + (j - i1)) % f1.length]) {
                match = false;
            }
        }
        if (match) {return true;}
    }
    return false;
}

Tuple!(ulong[], Apair) intersection_point(const Apair position, Apair direction, const Apair[] polygon) {
    real diam = 0;
    foreach (Apair p; polygon)
    {
        if ((position - p).abs() > diam) {
            diam = (position - p).abs();
        }
        foreach (Apair q; polygon)
        {
            if ((p - q).abs() > diam) {
                diam = (p - q).abs();
            }
        }
    }
    ulong[] index = [];
    auto dir = direction.unit() * diam;
    foreach (i,p ; polygon)
    {
        if (position.eq(p) || position.eq(polygon[(i+1) % polygon.length])) {
            if((position.eq(p) && polygon[(i+1) % polygon.length].eq(position + dir)) || 
                (polygon[(i+1) % polygon.length].eq(position) && p.eq(position + dir))) {
                    return tuple([i],Apair(0,0));
                }
            continue;
        }
        if (dist(p,polygon[(i+1) % polygon.length], position + dir) < 0.001) {
            index ~= i;
            continue;
        }
        auto intersects = lineSegmentLineSegmentIntersection(position, position + dir, p, polygon[(i+1) % polygon.length]);
        foreach(inter; intersects) {
            if ((position - inter).abs() < dir.abs()) {
                index = [i];
                dir = (inter - position);
            } else {
                index ~= i;
            }
        }
    }
    return tuple(index,dir);
}

real round(real x, uint places)
{
   real pwr = pow(10.0,places);
   return round(x * pwr) / pwr;
}

Apair round(Apair p) {
    return Apair(round(p.x,5), round(p.y,5));
}

real dist(Apair p1, Apair p2, Apair q) {
    assert(!p1.isnan());
    assert(!p2.isnan());
    return (Als(p1,p2).closestpoint(q) - q).abs();
}


extern (C++) {
    struct OffsetPolygon {
        double* coordinates;
        int length;
    }
    extern __gshared int buffer;

    OffsetPolygon offset(double p_Offset,ulong argc, double* argv , ulong hole_count, double* hole_points);
    struct Point {
        double x;
        double y;
        int id;
    }

    struct Edge {
        int start;
        int end;
    }

    struct StraightSkeleton {
        Point* points;
        int num_points;
        Edge* edges;
        int num_edges;
    }
    StraightSkeleton unweighted_skeleton(ulong outer_count, double *outer_points, ulong hole_count, double* hole_points);
    StraightSkeleton weighted_skeleton(ulong outer_count, double *outer_points);
}

Agraph!(2, Apair, void) generate_straight_skeleton(const Apair[] polygon_in, const Apair[] outer, ulong[] relevant_indices = [], 
        real factor = 1.0, bool simp = false, bool skeleton = false, bool ignore_simplification = false, bool weighted = true) {
    auto polygon = polygon_in.dup;
    if (!ignore_simplification) {
        polygon = simplify(polygon_in,outer, relevant_indices);
    }
    ulong[] relevant_edges = new ulong[](relevant_indices.length);
    relevant_edges[] = 0;
    bool nan = false;
    foreach(i,p ; polygon) {
        if (p.isnan()) {
            nan = true;
        }
        for(int j = 0; j < relevant_indices.length; j++){
            if (dist(p,polygon[(i+1) % polygon.length],polygon_in[relevant_indices[j]]) < 0.001 && !polygon[(i+1) % polygon.length].eq(polygon_in[relevant_indices[j]])) {
                relevant_edges[j] = i;
            }
        }
    }
    Apair[] relevant_positions = [];
    foreach(i ; relevant_indices) {
        relevant_positions ~= polygon_in[i].dup;
    }
    if (nan) {
        debug_output([tuple(polygon_in.dup,Ablue),tuple(outer.dup,Amagenta)],[], tuple(polygon,Ared));
        assert(false);
    }
    bool weight = (weighted && outer.length == 0 && relevant_indices.length == 2);
    double l1 = 1.0, l2 = 1.0;
    if (weight) {
        foreach(i,v ; polygon) {
            auto last = polygon[(i + polygon.length - 1) % polygon.length];
            if (i > relevant_edges[0] + 1 && i < relevant_edges[1]) {
                l1 += (v - last).abs();
            } else if (i > relevant_edges[0] + 1 && relevant_edges[0] > relevant_edges[1]) {
                l1 += (v - last).abs();
            } else if (i < relevant_edges[1] && relevant_edges[0] > relevant_edges[1]) {
                l1 += (v - last).abs();
            } else if(i == relevant_edges[0] + 1 && polygon[relevant_edges[0]].eq(polygon_in[relevant_indices[0]])) {
                l1 += (v - last).abs();
            } else if(i == relevant_edges[1] && polygon[relevant_edges[1]].eq(polygon_in[relevant_indices[1]])) {
                l1 += (v - last).abs();
            } else {
                l2 += (v - last).abs();
            }
        }
    }
    double[] polygon_args = new double[](buffer);
    auto diam = 0.0;
    foreach (i,v ; polygon) {
        polygon_args ~= v.x;
        polygon_args ~= v.y;
        if (weight) {
            if (i > relevant_edges[0] + 1 && i < relevant_edges[1]) {
                polygon_args ~= sqrt(sqrt(factor));
            } else if (i > relevant_edges[0] + 1 && relevant_edges[0] > relevant_edges[1]) {
                polygon_args ~= sqrt(sqrt(factor));
            } else if (i < relevant_edges[1] && relevant_edges[0] > relevant_edges[1]) {
                polygon_args ~= sqrt(sqrt(factor));
            } else if(i == relevant_edges[0] + 1 && polygon[relevant_edges[0]].eq(polygon_in[relevant_indices[0]])) {
                polygon_args ~= sqrt(sqrt(factor));
            } else if(i == relevant_edges[1] && polygon[relevant_edges[1]].eq(polygon_in[relevant_indices[1]])) {
                polygon_args ~= sqrt(sqrt(factor));
            } else {
                polygon_args ~= sqrt(sqrt(l2 / l1));
            }
        }
        foreach (w; polygon)
        {
            auto d = (v - w).abs();
            if (d > diam) {
                diam = d;
            }
        }
    }
    StraightSkeleton ss;
    if (outer.length > 0) {
        double[] border = new double[](buffer);
        for(ulong i = 0; i < outer.length; i++) {
            border ~= outer[i].x;
            border ~= outer[i].y;
        }
        ss = unweighted_skeleton(border.length, border.ptr, polygon_args.length, polygon_args.ptr);
    } else {
        if (weight) {
            ss = weighted_skeleton(polygon_args.length,polygon_args.ptr);
        } else {
            ss = unweighted_skeleton(polygon_args.length, polygon_args.ptr, 0, null);
        }
    }
    auto G = Agraph!(2,Apair,void)();
    for(int i = buffer;i < ss.num_points; i++) {
        auto pos = Apair(ss.points[i].x,ss.points[i].y);
        auto id = ss.points[i].id + 1;
        G.add(id, pos);
    }
    for(int i = buffer; i < ss.num_edges; i++) {
        auto s = ss.edges[i].start + 1;
        auto e = ss.edges[i].end + 1;
        if (s < e) {
            G.insert(s, e);
        }
    }
    if (skeleton) {
        Apath[] graph = [];
        foreach(u,v ; G) {
            graph ~= Als(G.val(u),G.val(v));
        }
        debug_output([tuple(polygon_in.dup,Ablue),tuple(polygon,Ared)], [tuple(graph,Amagenta)], tuple(relevant_positions,Ared));
    }
    // If we dont find start or endpoint, we assume it was removed because it is on a flat angle
    ulong[] filtered_relevant_indices = [];
    foreach(i ; relevant_indices) {
        if (!polygon.canFind(polygon_in[i])) {
            filtered_relevant_indices ~= i;
        }
    }
    foreach (ulong key; filtered_relevant_indices)
    {
        int index = G.vsup;
        G.add(index, polygon_in[key]);
        auto n = polygon_in[(key+1) % polygon_in.length];
        auto p = polygon_in[(key + polygon_in.length-1) % polygon_in.length];
        auto left = p - polygon_in[key];
        auto right = n - polygon_in[key];
        auto direction = bisector(left,right);
        if(weight) {
            auto angle = full_angle(right,left);
            auto w1 = sqrt(sqrt(factor));
            auto w2 = sqrt(sqrt(l2 / l1));
            auto prop = w1 / (w2 + w1);
            if (key == relevant_indices[0]) {
                prop = w2 / (w1 + w2);
            }
            direction = ARotate(angle * prop) * right.unit();
        }
        auto len = diam;
        Tuple!(int,int) edge = tuple(0,0);

        foreach(s,e ; G) {
            auto cur_intersection = lineSegmentLineSegmentIntersection(polygon_in[key],polygon_in[key] + len * direction, G.val(s), G.val(e));
            foreach(intersect ; cur_intersection) {
                edge = tuple(s,e);
                auto d = (intersect - polygon_in[key]).abs();
                if (d < len) {
                    len = d;
                }
            }
        }
        if (len == diam) {
            writeln(direction, polygon_in[key], key, polygon_in[key] + len * direction);
            break;
        }
        bool found = false;
        foreach (v ; G)
        {
            if (G.val(v).eq(polygon_in[key] + len * direction) && v != index) {
                found = true;
                G.insert(v,index);
            }
        }
        if (!found) {
            int n_index = G.vsup;
            G.add(n_index, polygon_in[key] + len * direction);
            G.insert(index,n_index);
            G.removeedge(edge[0],edge[1]);
            G.insert(edge[0],n_index);
            G.insert(n_index,edge[1]);
        }
    }
    return G;
}

Apair[] find_circle(const Apair[] polygon, real divisor = 3, bool simp = false, bool skeleton = false) {
    Apair[] relevant = simplify(polygon, [], []);
    double[] args = new double[](buffer);
    for(int i = 0;relevant.length > 2 && i < relevant.length; i++){
        args ~= to!double(relevant[i].x);
        args ~= to!double(relevant[i].y);
    }
    OffsetPolygon result = offset(divisor,args.length,args.ptr, 0, null);
    Apair[] result_poly = [];
    for (int i = buffer; i < result.length; i+=2) {
        auto p = Apair(result.coordinates[i],result.coordinates[i+1]);
        result_poly ~= round(p);
    }
    if (clockwise(relevant) != clockwise(result_poly)) {
        result_poly = result_poly.reverse;
    }
    return result_poly;
}

real[] place_along_path(Apath path, ulong amount) {
    real[] times = [0];
    for(int i = 1;i < (amount-1); i++) {
        real c = path.arctime(i * path.arclength / (amount - 1));
        times ~= c;
    }
    times ~= path.arctime(path.arclength);
    for(int i = 1; i + 1 < times.length; i++) {
        auto cur = path(round(times[i]));
        auto diff = (cur - path(times[i])).abs();
        if (diff < path.arclength / (2 * amount)) {
            times[i] = round(times[i]);
        }
    }
    return times;
}

Tuple!(Apair[],Apair[],Apair[],Apair[]) split_polygon(const Apair[] polygon, const Apair[] border) {
    auto p = polygon[0];
    auto intersection_data = first_intersection(0,polygon,border);
    auto index = intersection_data[0];
    auto found = intersection_data[1];
    auto intersect = p + intersection_data[2];
    auto target = p + 0.75 * (intersect - p);
    Apair[] polygon1, polygon2, outer1 = [], outer2 = [];
    if (!found) {
        polygon1 = [intersect, target] ~ polygon.dup ~ [polygon[0].dup,target, intersect];
        for(int i = 1; i <= border.length; i++) {
            if ((i != border.length && i != 1) || !(border[(index[0] + i) % border.length].eq(intersect,Aeps.e[3],0.001))){
                polygon1 ~= border[(index[0] + i) % border.length];
            }
        }
        polygon2 = [polygon[0].dup, target, intersect];
        for(int i = 1; i <= border.length; i++) {
            if ((i != border.length && i != 1) || !(border[(index[0] + i) % border.length].eq(intersect,Aeps.e[3],0.001))){
                polygon2 ~= border[(index[0] + i) % border.length];
            }
        }
        polygon2 ~= [intersect, target] ~ polygon.dup;
    } else {
        polygon1 = [intersect,target,polygon[0]];
        int last = 0;
        for(int i = 1; i < polygon.length; i++) {
            last = i;
            polygon1 ~= polygon[i];
            if(!polygon[i].eq(intersect,Aeps.e[3],0.001) && dist(polygon[i],polygon[(i+1) % polygon.length],intersect) < 0.001) {
                Apair next = polygon[(i + 2) % polygon.length];
                if (!polygon[(i+1) % polygon.length].eq(intersect,Aeps.e[3],0.001)) {
                    next = polygon[(i + 1) % polygon.length];
                }
                auto intersect_angle = full_angle(target - intersect,polygon[i] - intersect);
                auto next_angle = full_angle(next - intersect,polygon[i] - intersect);
                if (next_angle > intersect_angle) {
                    break;
                }
            } 
        }
        polygon2 = [polygon[0],target, intersect];
        for(ulong i = last + 1; i < polygon.length; i++) {
            polygon2 ~= polygon[i];
        }
    }
    if (polygon1.clockwise()) {
        outer1 = border.dup;
    }
    if (polygon2.clockwise()) {
        outer2 = border.dup;
    }
    return tuple(polygon1,outer1,polygon2,outer2);
}

Apair[] find_circle_with_articulation(const Apair[] polygon, const Apair[] border, bool simp = false, bool skeleton = false) {
    auto split = split_polygon(polygon,border);
    auto polygon1 = split[0];
    auto outer1 = split[1];
    auto polygon2 = split[2];
    auto outer2 = split[3];
    if(simp) {
        writeln("cycle  \n", polygon1,"\n",polygon2 ,"\n");
        writeln(simp,skeleton);
        debug_output([tuple(polygon.dup,Ablue),tuple(polygon1,Ared),tuple(polygon2,Agreen),tuple(border.dup,Amagenta)]);
    }
    auto path1 = find_polygon_path(polygon1, outer1, 2,1,1.0, simp, skeleton, false);
    auto path2 = find_polygon_path(polygon2, outer2, 1,0,1.0, simp, skeleton, false);
    auto total_path = path1 ~ path2[1..$];
    if (simp) {
        writeln("--- found --- ");
    }
    return total_path;
}

void draw_output(string name,Tuple!(Apair[],Apen)[] polygons, Tuple!(Apath[], Apen)[] drawings = [], Tuple!(Apair[],Apen) points = tuple(cast(Apair[])[], Ayellow)) {
    foreach(v ; polygons) {
        if (v[0].length == 0) {continue;}
        auto edges = Als(v[0][$-1],v[0][0]);
        for(int i = 0; i + 1 < v[0].length; i++) {
            edges ~= Als(v[0][i],v[0][i+1]);
        }
        draw(edges,v[1]);
        draw(v[0],v[1]);
    }
    foreach(v ; drawings) {
        foreach(p ; v[0]) {
            draw(p, v[1]);
        }
    }
    foreach(p ; points[0]) {
        draw(p,points[1]);
    }
    shiperase(name);
}

void debug_output(Tuple!(Apair[],Apen)[] polygons, Tuple!(Apath[], Apen)[] drawings = [], Tuple!(Apair[],Apen) points = tuple(cast(Apair[])[], Ayellow)) {
    foreach(v ; polygons) {
        if (v[0].length == 0) {continue;}
        auto edges = Als(v[0][$-1],v[0][0]);
        for(int i = 0; i + 1 < v[0].length; i++) {
            edges ~= Als(v[0][i],v[0][i+1]);
        }
        draw(edges,v[1]);
        draw(v[0],v[1]);
    }
    foreach(v ; drawings) {
        foreach(p ; v[0]) {
            draw(p, v[1]);
        }
    }
    draw(points[0],points[1]);
    shiperase("polygon_path_fail_total" ~ to!string(counter));
    foreach(v ; polygons) {
        if (v[0].length == 0) {continue;}
        auto edges = Als(v[0][$-1],v[0][0]);
        for(int i = 0;i + 1 < v[0].length; i++) {
            edges ~= Als(v[0][i],v[0][i+1]);
        }
        draw(edges,v[1]);
        draw(v[0],v[1]);
        draw(points[0],points[1]);
        shiperase("polygon_path_fail_partial" ~ to!string(++counter));
    }
}

real angle_factor = 4.0;

Tuple!(ulong[],bool,Apair) first_intersection(const int idx, const Apair[] polygon_in, const Apair[] border) {
    Apair[] polygon = [];
    for(int i = 0; i < polygon_in.length; i++) {
        auto cur = (idx + i) % polygon_in.length;
        if (polygon.length == 0 || !polygon_in[cur].eq(polygon[$-1],Aeps.e[3],0.001)) {
            if (i + 1 < polygon_in.length || !polygon_in[cur].eq(polygon[0],Aeps.e[3],0.001)) {
                polygon ~= polygon_in[cur];
            }
        }
    }
    auto prev = polygon[$-1];
    auto next = polygon[1];
    auto left = prev - polygon[0];
    auto right = next - polygon[0];
    auto bisect = bisector(left,right);
    auto direction = bisect;
    auto p = polygon[0];
    ulong[] index = [];
    bool polygon_intersection = false;
    auto angle = full_angle(right,left);
    real project = 0;
    real len = 0;
    for(int ind = 0; ind < polygon.length; ind++) {
        real cur_len;
        real cur_project;
        Apair cur_direction; 
        ulong[] cur_index = [];
        bool cur_polygon_intersect = true;
        if (polygon[ind].eq(p)) {
            if (ind > 0) {
                continue;
            }
            auto border_intersect = intersection_point(p,bisect,border);
            cur_len = border_intersect[1].abs();
            cur_project = cur_len;
            cur_direction = bisect;
            cur_index = border_intersect[0];
            cur_polygon_intersect = false;
        } else {
            cur_direction = (polygon[ind] - p).unit();
            if (cur_direction.eq(left.unit()) || cur_direction.eq(right.unit())) {
                continue;
            }
            cur_len = (polygon[ind] - p).abs();
            cur_project = (polygon[ind] - p) * bisect;
            if (cur_direction.eq(direction) && cur_project < project) {
                len = cur_len;
                project = cur_project;
            }
            cur_index = [ind, (ind + polygon.length - 1) % polygon.length];
        }
        auto dir_angle = full_angle(right,cur_direction);
        if (dir_angle <= angle / angle_factor || dir_angle >= angle * ((angle_factor - 1) / angle_factor)) {
            continue;
        }
        auto intersection = intersection_point(p,cur_direction,polygon);
        if (intersection[0].length > 0 && intersection[1].abs() <= cur_len) {
            cur_len = intersection[1].abs();
            cur_project = intersection[1] * bisect;
            cur_index = intersection[0];
            cur_polygon_intersect = true;
        }
        for(int i = 0;i < polygon.length; i++) {
            if (polygon[i].eq(p) || i == ind) {
                continue;
            }
            auto d = (polygon[i] - p).abs();
            auto a = full_angle(right,(polygon[i] - p).unit());
            auto projection = (polygon[i] - p) * cur_direction;
            auto lsquared = d * d - projection * projection;
            auto l = sqrt(max(lsquared,0));
            if (projection > 0.01 && d < cur_len && abs(a - dir_angle) < angle / 7.0 && l < 1.0) {
                cur_len = (polygon[i] - p).abs();
                cur_project = (polygon[i] - p) * bisect;
                cur_direction = (polygon[i] - p).unit();
                cur_index = [i, (i + polygon.length - 1) % polygon.length];
                cur_polygon_intersect = true;
            }
        }
        if (cur_project > project) {
            len = cur_len;
            project = cur_project;
            direction = cur_direction;
            index = cur_index;
            polygon_intersection = cur_polygon_intersect;
        }
    }
    if (len == 0) {
        writeln("result");
        project.writeln;
        polygon.writeln;
        direction.writeln;
        border.writeln;
    }
    assert(len > 0);
    return tuple(index,polygon_intersection,len * direction);
}

Apair[] simplify_remove_illegal(const Apair[] polygon, ulong[] indices) {
    if (polygon_area(polygon) < 0.001 || indices.length == 0) {
        return polygon.dup;
    }
    ulong index = indices[0];
    Apair[] intermediate = [];
    for(int i = 0; i < polygon.length; i++) {
        auto current = (i + indices[0]) % polygon.length;
        auto prev = (i + indices[0] + polygon.length - 1) % polygon.length;
        for(int j = 0; j < polygon.length; j++) {
            if (!polygon[j].eq(polygon[current],Aeps.e[3],0.001) && !polygon[j].eq(polygon[prev],Aeps.e[3],0.001) && dist(polygon[current],polygon[prev],polygon[j]) < 0.001){
                intermediate ~= polygon[j];
            }
        }
        if (i == 0) {
            index = intermediate.length;
        }
        if (intermediate.length == 0 || (!polygon[current].eq(intermediate[$-1],Aeps.e[3],0.001) && 
            (i + 1 < polygon.length || !polygon[current].eq(intermediate[0],Aeps.e[3],0.001)))){
            intermediate ~= polygon[current];
        }
    }
    Apair[] result = [];
    for(int i = 0; i < intermediate.length; i++) {
        auto current = (i + index) % intermediate.length;
        auto prev = (i + index + intermediate.length - 1) % intermediate.length;
        auto next = (i + index + 1) % intermediate.length;
        auto left = intermediate[prev] - intermediate[current];
        auto right = intermediate[next] - intermediate[current];
        auto angle = full_angle(right,left);
        result ~= intermediate[current];
        for(int j = i + 1; j < intermediate.length; j++) {
            auto cur = (j + index) % intermediate.length;
            if (current == cur || !intermediate[cur].eq(intermediate[current])) {continue;}
            auto r = intermediate[(cur + 1) % intermediate.length] - intermediate[current];
            if (full_angle(r,left) < angle) {
                i = j;
            }
        }
    }

    return result;
}

real simplify_factor = 2.0;

Apair move_by_line(int i, const Apair[] polygon,const Apair[] border) {
    auto prev = polygon[(i + polygon.length - 1) % polygon.length];
    auto next = polygon[(i + 1) % polygon.length];
    auto left = prev - polygon[i];
    auto right = next - polygon[i];
    auto bisector_intersect = prev + ((left.abs() / right.abs()) / (left.abs() / right.abs() + 1)) * (next - prev);
    auto bisector_len = (bisector_intersect - polygon[i]).abs();
    auto max_dist = bisector_len;
    auto bisect = bisector(left,right);
    if (left.unit().eq(-right.unit())) {
        max_dist = min(left.abs(),right.abs());
    }
    if (full_angle(right,left) >= 180) {
        foreach(j,p; polygon) {
            if (p.eq(polygon[i]) || polygon[(j + 1) % polygon.length].eq(polygon[i])) {
                continue;
            }
            auto intersect = lineSegmentLineSegmentIntersection(polygon[i],polygon[i] + max_dist * bisect, polygon[j],polygon[(j+1) % polygon.length]);
            foreach(inter ; intersect) {
                auto d = (inter - polygon[i]).abs();
                if (d < max_dist) {
                    max_dist = d;
                }
            }
        }
        foreach(j,p; border) {
            if (p.eq(polygon[i]) || border[(j + 1) % border.length].eq(polygon[i])) {
                continue;
            }
            auto intersect = lineSegmentLineSegmentIntersection(polygon[i],polygon[i] + max_dist * bisect, border[j],border[(j+1) % border.length]);
            foreach(inter ; intersect) {
                auto d = (inter - polygon[i]).abs();
                if (d < max_dist) {
                    max_dist = d;
                }
            }
        }
    }
    auto tri_max = max((polygon[i] - next).abs(),(polygon[i] - prev).abs(),(prev - next).abs());
    bool overlap = false;
    foreach (j,p; polygon)
    {
        if (p.eq(polygon[i]) || p.eq(prev) || p.eq(next) || !point_in_polygon([polygon[i],next,polygon[i] + max_dist * bisect,prev],p)) {
            if (i != j && p.eq(polygon[i])) {
                overlap = true;
            }
            continue;
        }
        auto intersect1 = lineSegmentLineSegmentIntersection(polygon[i],polygon[i] + max_dist * bisect,prev,prev + tri_max * (p - prev).unit());
        foreach(inter ; intersect1) {
            auto d = (inter - polygon[i]).abs();
            max_dist = min(d,max_dist);
        }
        auto intersect2 = lineSegmentLineSegmentIntersection(polygon[i],polygon[i] + max_dist * bisect,next,next + tri_max * (p - next).unit());
        foreach(inter ; intersect2) {
            auto d = (inter - polygon[i]).abs();
            max_dist = min(d,max_dist);
        }
    }
    if (!overlap) {
        return Apair(0,0);
    }
    foreach (p; border)
    {
        if (p.eq(polygon[i]) || p.eq(prev) || p.eq(next) || !point_in_polygon([polygon[i],next,polygon[i] + max_dist * bisect,prev],p)) {
            continue;
        }
        auto intersect1 = lineSegmentLineSegmentIntersection(polygon[i],polygon[i] + max_dist * bisect,prev,prev + tri_max * (p - prev).unit());
        foreach(inter ; intersect1) {
            auto d = (inter - polygon[i]).abs();
            max_dist = min(d,max_dist);
        }
        auto intersect2 = lineSegmentLineSegmentIntersection(polygon[i],polygon[i] + max_dist * bisect,next,next + tri_max * (p - next).unit());
        foreach(inter ; intersect2) {
            auto d = (inter - polygon[i]).abs();
            max_dist = min(d,max_dist);
        }
    }
    assert(max_dist > 0);
    return bisector(left,right) * max_dist / simplify_factor;
} 

Apair[] simplify_line(const Apair[] polygon,const Apair[] border, ulong[] indices) {
    Apair[] result = polygon.dup;
    for(int i = 0; i < result.length; i++) {
        if (indices.canFind(i)) {
            continue;
        }
        result[i] += move_by_line(i,result,border);   
    }
    return result;
}

Apair[2] move_by_angle(int i, const Apair[] polygon, const Apair[] border) {
    auto prev = polygon[(i + polygon.length - 1) % polygon.length];
    auto next = polygon[(i + 1) % polygon.length];
    auto left = prev - polygon[i];
    auto right = next - polygon[i];
    auto d1 = bisector(left,right).hat();
    auto d_left_len = left * d1;
    auto left_max = d_left_len * d1;
    auto d_right_len = right * d1;
    auto right_max = d_right_len * d1;
    foreach(j,p; polygon) {
        if (polygon[i].eq(p)) {
            continue;
        }
        if (!p.eq(prev)) {
            if (point_in_polygon([polygon[i],prev,polygon[i] + left_max], p)) {
                auto project = (p - prev) * (polygon[i] - prev).unit();
                auto l = sqrt(max((p - prev).abs() * (p - prev).abs() - project * project,0));
                auto L = l * (polygon[i] - prev).abs() / project;
                left_max = L * (d1);
            }
            
        }
        if(!polygon[i].eq(polygon[(j + 1) % polygon.length]) && 
            segmentsIntersect(polygon[j],polygon[(j + 1) % polygon.length], polygon[i], polygon[i] + left_max)) {
            auto intersect = lineSegmentLineSegmentIntersection(polygon[j],polygon[(j + 1) % polygon.length], polygon[i], polygon[i] + left_max);
            if ((intersect[0] - polygon[i]).abs() < left_max.abs()){
                left_max = (intersect[0] - polygon[i]);
            }
        }
        if (!p.eq(next)) {
            if (point_in_polygon([polygon[i],polygon[i] + right_max,next], p)) {
                auto project = (p - next) * (polygon[i] - next).unit();
                auto l = sqrt(max((p - next).abs() * (p - next).abs() - project * project,0));
                auto L = l * (polygon[i] - next).abs() / project;
                right_max = L * (-d1);
            }
        }
        if(!polygon[i].eq(polygon[(j + 1) % polygon.length]) && 
            segmentsIntersect(polygon[j],polygon[(j + 1) % polygon.length], polygon[i], polygon[i] + right_max)) {
            auto intersect = lineSegmentLineSegmentIntersection(polygon[j],polygon[(j + 1) % polygon.length], polygon[i], polygon[i] + right_max);
            if ((intersect[0] - polygon[i]).abs() < right_max.abs()){
                right_max = (intersect[0] - polygon[i]);
            }
        }
    }
    foreach(j,p; border) {
        if (point_in_polygon([polygon[i],prev,polygon[i] + left_max], p)) {
            auto project = (p - prev) * (polygon[i] - prev).unit();
            auto l = sqrt((p - prev).abs() * (p - prev).abs() - project * project);
            auto L = l * (polygon[i] - prev).abs() / project;
            left_max = L * (d1);
        }
        if(segmentsIntersect(border[j],border[(j + 1) % border.length], polygon[i], polygon[i] + left_max)) {
            auto intersect = lineSegmentLineSegmentIntersection(border[j],border[(j + 1) % border.length], polygon[i], polygon[i] + left_max)[0];
            left_max = (intersect - polygon[i]);
            if ((intersect - polygon[i]).abs() < left_max.abs()){
                left_max = (intersect - polygon[i]);
            }
        }
        if (point_in_polygon([polygon[i],polygon[i] + right_max,next], p)) {
            auto project = (p - next) * (polygon[i] - next).unit();
            auto l = sqrt((p - next).abs() * (p - next).abs() - project * project);
            auto L = l * (polygon[i] - next).abs() / project;
            right_max = L * (-d1);
        }
        if(segmentsIntersect(border[j],border[(j + 1) % border.length], polygon[i], polygon[i] + right_max)) {
            auto intersect = lineSegmentLineSegmentIntersection(border[j],border[(j + 1) % border.length], polygon[i], polygon[i] + right_max)[0];
            if ((intersect - polygon[i]).abs() < right_max.abs()){
                right_max = (intersect - polygon[i]);
            }
        }
    }
    assert(left_max * d1 > 0 && right_max * d1 < 0);
    return [left_max / simplify_factor,right_max / simplify_factor];
}

Apair[] simplify_angle(const Apair[] polygon, const Apair[] border) {
    Apair[] result = [];
    for(int i = 0; i < polygon.length; i++) {
        if (result.length > 0 && (result[$-1] - polygon[i]).abs() < 0.01) {
            continue;
        }
        auto left = polygon[(i + polygon.length - 1) % polygon.length] - polygon[i];
        auto right = polygon[(i+1) % polygon.length] - polygon[i];
        Apair p = polygon[i];
        auto a = full_angle(right,left);
        if (a > 270) {
            auto d = move_by_angle(i,polygon,border);
            result ~= p + d[0];
            result ~= p + d[1];
        } else if (abs(a - 180) > 2){
            result ~= p;
        }
    }
    if (polygon.length == 2 && !clockwise(result)) {
        result = result.reverse;
    }
    return result;
}

Apair[] simplify(const Apair[] polygon, const Apair[] border, ulong[] indices) {
    if (polygon.length <= 1) {
        return polygon.dup;
    }
    Apair[] removed = simplify_remove_illegal(polygon,indices.sort.array);
    ulong[] relevant_indices = new ulong[](indices.length);
    relevant_indices[] = removed.length;
    for(int i = 0; i < removed.length; i++) {
        for(int j = 0; j < indices.length; j++) {
            if (polygon[indices[j]].eq(removed[i])) {
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
    Apair[] simplified_lines = simplify_line(removed, border, relevant_indices);
    return simplify_angle(simplified_lines,border);
}
int counter = 0;


Apair[][] partition(Apair[] border, int[] sizes) {
    int remaining_total = 0;
    foreach (int v; sizes)
    {
        remaining_total += v;
    }
    Apair[][] result = [];
    if (sizes.length == 0) {
        return result;
    }
    Apair[] remaining = border.dup;
    while(sizes.length > 1) {
        auto factor = (sizes[$-1] + 0.0) / remaining_total;

        Apair dir = bisector(remaining[$-1] - remaining[0], remaining[1] - remaining[0]) * 500;
        for(int i = 1; i + 1 < remaining.length ; i++) {
            auto intersect = lineSegmentLineSegmentIntersection(remaining[0], remaining[0] + dir, remaining[i], remaining[i+1]);
            foreach(inter; intersect) {
                if ((remaining[0] - inter).abs() < dir.abs()) {
                    dir = inter - remaining[0];
                }
            }
        }
        int index = 0;
        foreach(i,v ; remaining) {
            foreach (j,w ; remaining) {
                auto d = (w - v);
                if (((abs(to!int(i) - to!int(j)) > 1) && (i * j > 0 || (i + j + 1) < remaining.length)) && d.abs() > dir.abs()) {
                    index = to!int(i);
                    dir = d;
                }
            }
        }
        auto line = dir * factor;
        Apair p = remaining[index] + line;
        auto w = ARotate(90) * dir;
        auto l1 = w;
        auto l2 = -w;
        auto index1 = -1;
        auto index2 = -1;
        for(int i = 0; i < remaining.length; i++) {
            auto intersections1 = lineSegmentLineSegmentIntersection(p,p + l1, remaining[i], remaining[(i+1) % remaining.length]);
            auto intersections2 = lineSegmentLineSegmentIntersection(p,p + l2, remaining[i], remaining[(i+1) % remaining.length]);
            foreach (Apair intersect; intersections1)
            {
                auto l = (intersect - p);
                index1 = i;
                if (l.abs() < l1.abs()) {
                    l1 = l;
                }
            }
            foreach (Apair intersect; intersections2)
            {
                auto l = (intersect - p);
                index2 = i;
                if (l.abs() < l2.abs()) {
                    l2 = l;
                }
            }
        }
        Apair[] part = [p + l1];
        ulong current = (index1 + 1) % remaining.length;
        while(current != (index2 + 1) % remaining.length) {
            if (!part[$-1].eq(remaining[current])){
                part ~= remaining[current];
            }
            current = (current + 1) % remaining.length;
        }
        if (!part[$-1].eq(p + l2)){
            part ~= p + l2;
        }
        Apair[] rem = [p + l2];
        while(current != (index1 + 1) % remaining.length) {
            if (!rem[$-1].eq(remaining[current])){
                rem ~= remaining[current];
            }
            current = (current + 1) % remaining.length;
        }
        if (!rem[$-1].eq(p + l1)){
            rem ~= p + l1;
        }
        result = part ~ result;
        remaining = rem;
        remaining_total -= sizes[$-1];
        sizes = sizes[0..$-1];
    }
    result = remaining ~ result;
    return result;
}
import std.random : Random, uniform;
auto rnd = Random(42);

Apair[] find_polygon_path(const Apair[] polygon_in, const Apair[] outer, ulong start_index_in, ulong end_index_in, 
        real factor = 1.0, bool simp = false, bool skeleton = false, bool weighted = true) {
    auto G = generate_straight_skeleton(polygon_in,outer,[start_index_in,end_index_in],factor,simp,skeleton,false, weighted);
    auto start = 0,end = 0;
    auto min_start = 0.001;
    auto min_end = 0.001;
    foreach(v ; G) {
        auto p = G.val(v);
        if (p.eq(polygon_in[start_index_in],Aeps.e[3],min_start)) {
            start = v;
            min_start = (p - polygon_in[start_index_in]).abs();
        }
        if (p.eq(polygon_in[end_index_in],Aeps.e[3],min_end)) {
            end = v;
            min_end = (p - polygon_in[end_index_in]).abs();
        }
    }
    Apair[] path = [];
    foreach(v;G.Avminpath(start,end)) {
        auto p = G.val(v);
        if (path.length == 0 || !p.eq(path[$-1])){
            path ~= G.val(v);
        }
    }
    auto fail = path.length == 0;
    for(int i = 1; !fail && i+1 < path.length; i++) {
        for(int j = 0; j < polygon_in.length; j++) {
            if (j == start_index_in || j == end_index_in) {continue;}
            if (path[i].eq(polygon_in[j])) {fail = true;}
        }
        for(int j = 0; j < outer.length; j++) {
            if (path[i].eq(outer[j])) {fail = true;}
        }
    }
    if (fail || path.length < 2 || skeleton) {
        Apath[] path_poly = [];
        for(int i = 0;i + 1 < path.length; i++) {
            path_poly ~= Als(path[i],path[i+1]);
        }
        Apath[] graph = [];
        foreach(u,v ; G) {
            graph ~= Als(G.val(u),G.val(v));
        }
        writeln("outputting skeleton at ",counter);
        writeln(start_index_in, "," ,end_index_in);
        Apair[] polygon = simplify(polygon_in,outer, [start_index_in,end_index_in]);
        debug_output([tuple(polygon_in.dup,Ayellow),tuple(polygon.dup,Ablue), tuple(outer.dup, Ablue)], [tuple(graph,Amagenta),tuple(path_poly,Agreen)]);
        if (fail || path.length < 2) {
            writeln("ERROR finding legal path in generated straight skeleton");
            polygon_in.writeln;
            writeln(polygon_in[start_index_in],polygon_in[end_index_in]);
            outer.writeln;
            G.writeln;
            path.writeln;
            assert(false);
        }
    }
    //polygon.writeln;

    Apair[] result = [polygon_in[start_index_in]];
    for (int i = 1;i + 1 < path.length; i++) {
        if (path[i].eq(polygon_in[end_index_in],Aeps.e[3],0.001)) {
            break;
        }
        if (!path[i].eq(result[$-1],Aeps.e[3],0.001)) {
            result ~= round(path[i]);
        }
    }
    if (!result[$-1].eq(polygon_in[end_index_in],Aeps.e[3],0.001)){
        result ~= polygon_in[end_index_in];
    }
    return result;
}

Apair bisector(const Apair left, const Apair right) {
    assert(!right.iszero() && !left.iszero());
    auto r = right.unit();
    auto l = left.unit();
    auto bisector = l + r;
    if (l.eq(r) || full_angle(r,l) > 180.0) {
        bisector = -bisector;
    }else if(bisector.iszero()) {
        bisector = ARotate(90) * r;
    }
    return bisector.unit();
}

Apair[] outer_circle(Apair[] original) {
    auto maxX = original[0].x, minX = original[0].x;
    auto maxY = original[0].y, minY = original[0].y;
    auto sum = Apair(0,0);
    for(int i = 1; i < original.length; i++) {
        if (original[i].x > maxX) {
            maxX = original[i].x;
        }
        if (original[i].x < minX) {
            minX = original[i].x;
        }
        if (original[i].y > maxY) {
            maxY = original[i].y;
        }
        if (original[i].y < minY) {
            minY = original[i].y;
        }
        sum += original[i];
    }

    auto center = sum / original.length;
    Apair[] result = [];
    auto d = 0.0;
    for(int i = 0; i < original.length; i++) {
        if ((original[i] - center).abs() > d) {
            d = (original[i] - center).abs();
        } 
    }
    auto l = (d * 2 + original.length * size);
    auto s = original.length;
    for(ulong i = 0; i < s; i++) {
        result ~= center + ARotate(i * 360.0 / s) * Apair(0,l);
    }
    return result; 
}

bool clockwise(const Apair[] polygon) {
    auto sum = 0.0;
    for(int i = 0; i < polygon.length; i++) {
        Apair p0 = polygon[i];
        Apair p1 = polygon[(i+1)%polygon.length];
        sum += (p1.x - p0.x) * (p1.y + p0.y);
    }
    return sum > 0;
}

real polygon_area(const Apair[] polygon) {
    auto poly = clockwise(polygon) ? polygon.dup.reverse : polygon;
    real result = 0;
    for(int i = 0; i < poly.length; i++) {
        result += poly[i].x * (poly[(i+1) % poly.length].y - poly[(i-1 + poly.length) % poly.length].y);
    }
    return result * 0.5;
}

Apair[] intersection_paths(Apath p1, Apath p2) {
    Apair[] intersections = [];
    for(int i = 0; i < p1.arctime(p1.arclength); i++) {
        for(int j = 0; j < p2.arctime(p2.arclength); j++) {
            auto intersect = lineSegmentLineSegmentIntersection(round(p1(i)), round(p1(i+1)), round(p2(j)), round(p2(j+1)));
            //if (intersect.length > 0) {
            //    intersect.writeln;
            //    writeln(to!string(p1(i)) ~ "." ~ to!string(p1(i+1)) ~ "--" ~ to!string(p2(j)) ~ "." ~ to!string(p2(j+1)) ~ "\n");
            //}
            intersections ~= intersect;
            
        }
    }
    return intersections;
}

real full_angle(const Apair p, const Apair q) {
    assert(!p.iszero() && !q.iszero());
    if (p.unit().eq(q.unit())) {
        return 360;
    }
    if (p.unit().eq(-q.unit())) {
        return 180;
    }
    real dot = p.x * q.x + p.y * q.y;
    real det = p.x * q.y - p.y * q.x;
    real result = atan2(det,dot);
    if (result < 0) {
        result += 2 * PI;
    }
    return (result/(2 * PI)) * 360;
}

bool point_in_polygon(const Apair[] polygon, const Apair point) {
    ulong num_vertices = polygon.length;
    real x = point.x, y = point.y;
    bool inside = false;

    // Store the first point in the polygon and initialize
    // the second point
    Apair p1 = polygon[0], p2;

    // Loop through each edge in the polygon
    for (int i = 1; i <= num_vertices; i++) {
        // Get the next point in the polygon
        p2 = polygon[i % num_vertices];
        if(p2.eq(point)) {
            i++;
            p1 = polygon[i % num_vertices];
            continue;
        }
        // Check if the point is above the minimum y
        // coordinate of the edge
        if (y > min(p1.y, p2.y)) {
            // Check if the point is below the maximum y
            // coordinate of the edge
            if (y <= max(p1.y, p2.y)) {
                // Check if the point is to the left of the
                // maximum x coordinate of the edge
                // Calculate the x-intersection of the
                // line connecting the point to the edge
                if (x <= max(p1.x, p2.x)) {
                    // Check if the point is on the same
                    // line as the edge or to the left of
                    // the x-intersection
                    real x_intersection
                        = (y - p1.y) * (p2.x - p1.x)
                            / (p2.y - p1.y)
                        + p1.x;
                    if (p1.x == p2.x
                        || x <= x_intersection) {
                        // Flip the inside flag
                        inside = !inside;
                    }
                }
            }
        }

        // Store the current point as the first point for
        // the next iteration
        p1 = p2;
    }

    // Return the value of the inside flag
    return inside;
}