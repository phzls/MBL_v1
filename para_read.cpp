//
// Created by Liangsheng Zhang on 4/17/15.
//

#include <fstream>
#include "parameters.h"
#include "para_read.h"

using namespace std;

void para_get(string filename, const vector<vector<string> >& content, string keyword, int& var){
    string temp;
    bool find = false;
    for (int i=0;i<content.size();i++){
        if (content[i][0] == keyword){
            find = true;
            temp = content[i][1];
            break;
        }
    }

    if (!find){
        cout << keyword << " in " << filename << ".dat cannot be found." << endl;
        abort();
    }

    stringstream convert(temp);
    convert >> var;
}

void para_get(string filename, const vector<vector<string> >& content, string keyword, double& var){
    string temp;
    bool find = false;
    for (int i=0;i<content.size();i++){
        if (content[i][0] == keyword){
            find = true;
            temp = content[i][1];
            break;
        }
    }

    if (!find){
        cout << keyword << " in " << filename << ".dat cannot be found." << endl;
        abort();
    }

    stringstream convert(temp);
    convert >> var;
}

void para_get(string filename, const vector<vector<string> >& content, string keyword, string& var){
    bool find = false;
    for (int i=0;i<content.size();i++){
        if (content[i][0] == keyword){
            find = true;
            var = content[i][1];
            break;
        }
    }

    if (!find){
        cout << keyword << " in " << filename << ".dat cannot be found." << endl;
        abort();
    }
}

void para_get(string filename, const vector<vector<string> >& content, string keyword, bool& var){
    string temp;
    bool find = false;
    for (int i=0;i<content.size();i++){
        if (content[i][0] == keyword){
            find = true;
            temp = content[i][1];
            break;
        }
    }

    if (!find){
        cout << keyword << " in " << filename << ".dat cannot be found." << endl;
        abort();
    }

    if (temp == "true") var = true;
    else if (temp == "false") var = false;
    else{
        cout << keyword << " in " << filename << ".dat is expected to be a bool." << endl;
        cout << keyword << " value: " << temp << endl;
        abort();
    }
}

void para_file_read(string filename, vector<vector<string> >& content){

    stringstream file;
    file << "./parameters/" << filename << ".dat";

    ifstream fin(file.str().c_str());

    if (!fin.good()){
        cout << filename << ".dat is not accessible." << endl;
        abort();
    }

    int index = 0;
    while(!fin.eof()){
        string temp_line;
        getline(fin,temp_line);

        if (temp_line.substr(0,2) != "//" && temp_line.size() >1){
            // If it is not a comment line

            int c_start = temp_line.find("//"); // The starting positon of the comment

            string short_line = temp_line;
            if (c_start != string::npos){
                // comment exists
                short_line = temp_line.substr(0,c_start);
            }

            istringstream buf(short_line);
            istream_iterator<string> beg(buf), end;
            vector<string> tokens(beg, end); // Split into works by blank space

            if (tokens.size() != 2){
                cout << "The line in " << filename << ".dat has an invalid line whose number of"
                     << " variable is not 2." << endl;
                cout << "Line: " << short_line << endl;
                abort();
            }

            content.push_back(tokens);
        }
    }
}
