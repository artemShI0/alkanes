    for(int i = 1; i < 17; ++i){
        string s = ".\\skelets\\C";
        s += char(i / 10 + '0');
        s += char(i % 10 + '0');
        CreateDirectoryA(s.c_str(), NULL);
    }