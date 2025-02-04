                    graph[graph[v][i]].x = graph[v].x + L * cos((-2 * alpha + fi) * (!st) + (2 * alpha - _fi) * (st));
                    graph[graph[v][i]].y = graph[v].y + L * sin((-2 * alpha + fi) * (!st) + (2 * alpha - _fi) * (st));
                    coordinate_dfs(graph[v][i], v, !st, diameter, dist + 1);