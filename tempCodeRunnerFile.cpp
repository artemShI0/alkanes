        // for (int j = 0; j < molecules[i].size(); j += m){
        //     vector<thread> th;
        //     for(int l = 0; l < m && j + l < molecules[i].size(); ++l){
        //         thread t(&Molecule::build, molecules[i][j + m]);
        //         th.push_back(move(t)); 
        //     }
        //     for(int l = 0; l < m && j + l < molecules[i - 1].size(); ++l){
        //         th[l].join();
        //     }  
        // }