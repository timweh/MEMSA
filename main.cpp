#include <iostream>
#include <fstream>
#include <chrono>
using namespace std;

string directory;
string ref_file_name;
string read_file_name;
string res_file_name;

int min_seed_size = 20;
int max_merge_gap = 1;

class Seed {
    private:
        int reference_index;
        int read_index;
        int length;

    public:
        Seed(string indices) {
            int delim_pos = indices.find("\t");
            reference_index = stoi(indices.substr(0, delim_pos))-1;
            indices = indices.substr(delim_pos+1);

            delim_pos = indices.find("\t");
            read_index = stoi(indices.substr(0, delim_pos))-1;
            indices = indices.substr(delim_pos+1);

            length = stoi(indices);
        }

        Seed(int ref, int read, int len) {
            reference_index = ref;
            read_index = read;
            length = len;
        }

        int ref() {
            return reference_index;
        }

        int read() {
            return read_index;
        }

        int len() {
            return length;
        }

        int ref_end() {
            return reference_index + length - 1;
        }

        int read_end() {
            return read_index + length - 1;
        }

        bool contains(int refIndex) {
            return ((refIndex <= ref_end()) && (refIndex >= reference_index));
        }

        void setStart(int refIndex) {
            length += reference_index - refIndex;
            read_index += refIndex - reference_index;
            reference_index = refIndex;
        }

        void setEnd(int refIndex) {
            length += refIndex - reference_index;
        }

        void cutStart(int cutOffCharacters) {
            read_index += cutOffCharacters;
            reference_index += cutOffCharacters;
            length -= cutOffCharacters;
        }

        void cutEnd(int cutOffCharacters) {
            length -= cutOffCharacters;
        }
};

class ReadMEMs {
    private:
        string readName;
        vector<Seed> seeds;

    public:
        ReadMEMs(string name) {
            readName = name;
        }

        void add_seed(Seed seed) {
            if(seeds.size() > 0) {
                int remove = 0;
                if(seed.read_end() >= seeds.back().read()) {
                    remove = seed.read_end() - seeds.back().read() + 1;
                    seed.cutEnd(remove);
                    seeds.back().cutStart(remove);
                }
                if(seed.ref_end() >= seeds.back().ref()) {
                    remove = seed.ref_end() - seeds.back().ref() + 1;
                    seed.cutEnd(remove);
                    seeds.back().cutStart(remove);
                }
            }
            seeds.push_back(seed);
        }

        int seed_count() {
            return seeds.size();
        }

        Seed get_seed(int index) {
            return seeds.at(index);
        }

        int ref_index(int read_index) {
            int i = 0;
            while(read_index < seeds.at(i).read()) {
                i++;
            }
            return read_index + seeds.at(i).ref() - seeds.at(i).read();
        }

        int read_index(int ref_index) {
            int i = 0;
            while(ref_index < seeds.at(i).ref()) {
                i++;
            }
            if(ref_index <= seeds.at(i).ref_end()) {
                return ref_index + seeds.at(i).read() - seeds.at(i).ref();
            } else {
                cerr << "[ERROR] INDEX CONVERSION FAILED - REF_INDEX IS NOT CONTAINED IN ANY SEED\n" << endl;
                throw 1;
            }
        }

        void intersect(vector<int>& start, vector<int>& end) {
            vector<int> new_start;
            vector<int> new_end;
            int thisIndex = 0;
            int listIndex = 0;
            while (thisIndex < seeds.size() && listIndex < start.size()) {
                if (!(start.at(listIndex) > seeds.at(thisIndex).ref_end() || seeds.at(thisIndex).ref() > end.at(listIndex))) {
                    new_start.push_back(max(start.at(listIndex), seeds.at(thisIndex).ref()));
                    new_end.push_back(min(end.at(listIndex), seeds.at(thisIndex).ref_end()));
                }
                if (start.at(listIndex) > seeds.at(thisIndex).ref()) {
                    listIndex++;
                } else {
                    thisIndex++;
                }
            }
            start = new_start;
            end = new_end;
        }

        void clearSeeds() {
            seeds.clear();
        }
};

void readFasta(string file_name, vector<string>& seq_names, vector<string>& sequences) {
    fstream fileInput;
    fileInput.open(file_name, ios::in);
    if (fileInput.is_open()) {
        string line;
        while (getline(fileInput, line)) {
            if (!line.empty() && line.at(0) == '>') {
                seq_names.push_back(line);
            } else if (!line.empty()) {
                if (seq_names.size() > sequences.size()) {
                    sequences.push_back(line);
                } else {
                    sequences.back().append(line);
                }
            }
        }
        fileInput.close();
    } else {
        cerr << "[ERROR] FILE COULD NOT BE OPENED\n" << endl;
        throw 1;
    }
}

void readFasta(string file_name, vector<string>& sequences) {
    vector<string> names;
    readFasta(file_name, names, sequences);
}

void writeFasta(string file_name, vector<string>& seq_names, vector<string>& sequences) {
    fstream fileOutput;
    fileOutput.open(file_name, ios::out);
    if (fileOutput.is_open()) {
        for (int i = 0; i < sequences.size(); i++) {
            fileOutput << seq_names.at(i) << endl;
            fileOutput << sequences.at(i) << endl;
        }
        fileOutput.close();
    } else {
        cerr << "[ERROR] FILE COULD NOT BE OPENED\n" << endl;
        throw 1;
    }
}

vector<ReadMEMs> reads;
vector<string> sequences;

void extractSubsequences(int ref_start, int ref_end, vector<string>& subseq, bool seed) {
    if (ref_end < ref_start && ref_end != -1) {
        cerr << "[ERROR] END INDEX IS SMALLER THAN START\n" << endl;
        throw 1;
    }
    subseq.clear();
    int start_index = ref_start;
    int end_index = ref_end;
    if (!seed) {
        start_index++;
        end_index--;
    }
    if (ref_start == -1) {
        start_index = 0;
    }
    if (ref_end == -1) {
        end_index = sequences.at(0).size()-1;
    }
    subseq.push_back(sequences.at(0).substr(start_index, end_index-start_index+1));
    for (int i = 0; i < reads.size(); i++) {
        if (seed) {
            start_index = reads.at(i).read_index(ref_start);
            end_index = reads.at(i).read_index(ref_end);
        } else {
            if (ref_start == -1) {
                start_index = 0;
            } else {
                start_index = reads.at(i).read_index(ref_start)+1;
            }
            if (ref_end == -1) {
                end_index = sequences.at(i+1).size()-1;
            } else {
                end_index = reads.at(i).read_index(ref_end)-1;
            }
        }
        subseq.push_back(sequences.at(i+1).substr(start_index, end_index-start_index+1));
        if (end_index-start_index+1 < 0) {
            cerr << "[ERROR] CANNOT EXTRACT NEGATIVE SUBSEQUENCE\n" << endl;
            throw 1;
        }
    }
}

void alignmentScore(vector<string>& alignment, int& gaps, int& matches, int& mismatches, int& col_matches) {
    gaps = 0;
    matches = 0;
    mismatches = 0;
    col_matches = 0;
    bool hasGap;
    bool hasMismatch;
    int a,c,g,t;
    for (int i = 0; i < alignment.at(0).length(); i++) {
        hasGap = false;
        hasMismatch = false;
        a = 0;
        c = 0;
        g = 0;
        t = 0;
        for (int j = 0; j < alignment.size(); j++) {
            switch (toupper(alignment.at(j).at(i))) {
                case 'A': a++; break;
                case 'C': c++; break;
                case 'G': g++; break;
                case 'T': t++; break;
                case 'M': a++; c++; break;
                case 'R': a++; g++; break;
                case 'W': a++; t++; break;
                case 'S': c++; g++; break;
                case 'Y': c++; t++; break;
                case 'K': g++; t++; break;
                case 'V': a++; c++; g++; break;
                case 'H': a++; c++; t++; break;
                case 'D': a++; g++; t++; break;
                case 'B': c++; g++; t++; break;
                case 'N': a++; c++; t++; g++; break;
                default: break;
            }
            if (alignment.at(j).at(i) != alignment.at(0).at(i)) {
                hasMismatch = true;
            }
            if (alignment.at(j).at(i) == '-') {
                hasGap = true;
                gaps++;
            }
        }
        col_matches += max(max(max(a,c),t),g)-1;
        if (hasMismatch && !hasGap) {
            mismatches++;
        } else if (!hasGap) {
            matches++;
        }
    }
}

void writeResult(string file_name, int exec_time) {
    fstream fileOutput;
    fileOutput.open("results.csv", ios::app);
    if (fileOutput.is_open()) {
        fileOutput << file_name << "," << min_seed_size << "," << max_merge_gap << "," << exec_time << endl;
        fileOutput.close();
    } else {
        cerr << "[ERROR] FILE COULD NOT BE OPENED\n" << endl;
        throw 1;
    }
}

void writeAnalysis(string file_name, double mmcr, double ampc) {
    fstream fileOutput;
    fileOutput.open("analysis.csv", ios::app);
    if (fileOutput.is_open()) {
        fileOutput << file_name << "," << mmcr << "," << ampc << endl;
        fileOutput.close();
    } else {
        cerr << "[ERROR] FILE COULD NOT BE OPENED\n" << endl;
        throw 1;
    }
}

void cutLength(string filename, string file_ending, int length) {
    vector<string> names;
    vector<string> reads;
    readFasta(filename+"."+file_ending, names, reads);
    while(names.size() >= length) {
        names.pop_back();
        reads.pop_back();
    }
    writeFasta(filename+"-"+to_string(length)+"."+file_ending, names, reads);
}

void analyzeResults(string filename) {
    vector<string> alignment;
    int gaps, matches, mismatches, col_matches;

    readFasta(filename, alignment);
    alignmentScore(alignment, gaps, matches, mismatches, col_matches);

    cout << filename << endl;

    double match_ratio = matches;
    match_ratio /= (mismatches+matches);
    double matches_per_col = col_matches;
    matches_per_col /= alignment.at(0).length();
    matches_per_col /= alignment.size()-1;

    cout << "MMCR: " << match_ratio << endl;
    cout << "AMCP: " << matches_per_col << endl;

    writeAnalysis(filename,match_ratio,matches_per_col);
}

void performMafft() {
    // ---------- CALL MSA FOR COMPARISON ----------

    // --- COMBINE INPUT FILES INTO SINGLE FILE ---

    ifstream if_ref(ref_file_name, ios_base::binary);
    ifstream if_read(read_file_name, ios_base::binary);
    ofstream of_all(directory + "temp/all.fa", ios_base::binary);
    of_all << if_ref.rdbuf() << endl << if_read.rdbuf();
    if_ref.close();
    if_read.close();
    of_all.close();

    // --- EXECUTE MAFFT ---
    system((directory + "mafft/mafft.bat " + directory + "temp/all.fa > " + res_file_name).c_str());
}

void performMEMSA() {
    // ---------- CALL slaMEM TO GENERATE MEM SEEDS ----------

    system((directory + "slaMEM/slaMEM -l " + to_string(min_seed_size) + " -o " + directory + "temp/temp_mems.txt " + ref_file_name + " " + read_file_name).c_str());

    // ---------- LOAD MEM INDICES FROM TXT FILE ----------

    fstream indexFileInput;
    indexFileInput.open(directory + "temp/temp_mems.txt", ios::in);
    if(indexFileInput.is_open()) {
        string line;
        while (getline(indexFileInput, line)) {
            if (line.at(0) == '>') {
                reads.push_back(ReadMEMs(line));
            } else if (!line.empty()) {
                reads.back().add_seed(Seed(line));
            }
        }
        indexFileInput.close();
    }

    // ---------- FIND COMMON SEEDS ----------

    vector<int> start;
    vector<int> end;

    for(int i = 0; i < reads.at(0).seed_count(); i++) {
        start.push_back(reads.at(0).get_seed(i).ref());
        end.push_back(reads.at(0).get_seed(i).ref_end());
    }

    int total_seeds = reads.at(0).seed_count();
    for(int i = 1; i < reads.size(); i++) {
        reads.at(i).intersect(start, end);
        total_seeds += reads.at(i).seed_count();
    }

    if (start.empty()) {
        cerr << "[ERROR] NO COMMON SEEDS FOUND" << endl;
        throw 0;
    }

    // ---------- LOAD SEQUENCES FROM FASTA FILES ----------

    vector<string> seq_names;
    readFasta(ref_file_name, seq_names, sequences);
    readFasta(read_file_name, seq_names, sequences);

    // ---------- MERGE SEEDS WITH GAP BELOW THRESHOLD ----------

    int ref_gap;
    int read_gap;
    bool equal_gap;

    vector<int> mstart;
    vector<int> mend;

    // --- MERGE LAST SEED WITH END ---
    ref_gap = sequences.at(0).length() - end.at(0) - 1;
    equal_gap = true;
    if (ref_gap <= max_merge_gap) {
        for(int i = 0; i < reads.size(); i++) {
            read_gap = sequences.at(i+1).length() - reads.at(i).read_index(end.at(0));
            if (read_gap != ref_gap) {
                equal_gap = false;
                break;
            }
        }
        if(equal_gap) {
            mend.push_back(sequences.at(0).length()-1);
        } else {
            mend.push_back(end.at(0));
        }
    } else { // merge with end
        mend.push_back(end.at(0));
    }
    mstart.push_back(start.at(0));

    // --- MERGE SEEDS WITH EACH OTHER
    for (int i = 1; i < start.size(); i++) {
        ref_gap = start.at(i-1) - end.at(i) - 1;
        equal_gap = true;
        if (ref_gap <= max_merge_gap) {
            for(int j = 0; j < reads.size(); j++) {
                read_gap = reads.at(j).read_index(start.at(i-1))-reads.at(j).read_index(end.at(i))-1;
                if (read_gap != ref_gap) {
                    equal_gap = false;
                    break;
                }
            }
            if(equal_gap) {
                mstart.back() = start.at(i);
            } else {
                mstart.push_back(start.at(i));
                mend.push_back(end.at(i));
            }
        } else {
            mstart.push_back(start.at(i));
            mend.push_back(end.at(i));
        }
    }

    // --- MERGE FIRST SEED WITH START ---
    ref_gap = mstart.back();
    equal_gap = true;
    if (ref_gap <= max_merge_gap) { // merge with start
        for(int i = 0; i < reads.size(); i++) {
            read_gap = reads.at(i).read_index(mstart.back());
            if (read_gap != ref_gap) {
                equal_gap = false;
                break;
            }
        }
        if(equal_gap) {
            mstart.back() = 0;
        }
    }

    cout << endl << total_seeds << " -> " << start.size() << " -> " << mstart.size() << " seeds" << endl << endl;

    // ---------- PERFORM MSA ON SUBSEQUENCES BETWEEN SEEDS ----------

    vector<string> alignments(reads.size()+1);
    vector<string> sseq;

    int start_index;
    int end_index;

    // --- EXTRACT PRE-SEED SEQUENCE ---
    if (mstart.back() == 0) {
        for (int i = 0; i < reads.size(); i++) {
            if (reads.at(i).read_index(mstart.back()) != 0) {
                extractSubsequences(-1, mstart.back(), sseq, false);
                // --- CREATE SUBSEQUENCE FASTA FILES ---
                writeFasta(directory + "temp/temp_sseq.fa", seq_names, sseq);

                // --- CALL MSA ALGORITHM ---
                system((directory + "mafft/mafft.bat " + directory + "temp/temp_sseq.fa > " + directory + "temp/temp_msa.fa").c_str());

                // --- READ MSA OUTPUT FILE ---
                sseq.clear();
                readFasta(directory + "temp/temp_msa.fa", sseq);
                for(int j = 0; j < sseq.size(); j++) {
                    transform(sseq.at(j).begin(), sseq.at(j).end(), sseq.at(j).begin(), ::toupper);
                    alignments.at(j).append(sseq.at(j));
                }
                break;
            }
        }
    } else {
        extractSubsequences(-1, mstart.back(), sseq, false);
        // --- CREATE SUBSEQUENCE FASTA FILES ---
        writeFasta(directory + "temp/temp_sseq.fa", seq_names, sseq);

        // --- CALL MSA ALGORITHM ---
        system((directory + "mafft/mafft.bat " + directory + "temp/temp_sseq.fa > " + directory + "temp/temp_msa.fa").c_str());

        // --- READ MSA OUTPUT FILE ---
        sseq.clear();
        readFasta(directory + "temp/temp_msa.fa", sseq);
        for(int j = 0; j < sseq.size(); j++) {
            transform(sseq.at(j).begin(), sseq.at(j).end(), sseq.at(j).begin(), ::toupper);
            alignments.at(j).append(sseq.at(j));
        }
    }

    for (int i = mstart.size()-1; i >= 0; i--) {

        // --- EXTRACT SEEDS ---
        extractSubsequences(mstart.at(i), mend.at(i), sseq, true);
        for(int j = 0; j < sseq.size(); j++) {
            alignments.at(j).append(sseq.at(j));
        }

        if (i == 0) {
            break;
        }

        // --- EXTRACT SUBSEQUENCES ---
        extractSubsequences(mend.at(i), mstart.at(i-1), sseq, false);

        // --- CREATE SUBSEQUENCE FASTA FILES ---
        writeFasta(directory + "temp/temp_sseq.fa", seq_names, sseq);

        // --- CALL MSA ALGORITHM ---
        system((directory + "mafft/mafft.bat " + directory + "temp/temp_sseq.fa > " + directory + "temp/temp_msa.fa").c_str());

        // --- READ MSA OUTPUT FILE ---
        sseq.clear();
        readFasta("temp/temp_msa.fa", sseq);
        for(int j = 0; j < sseq.size(); j++) {
            transform(sseq.at(j).begin(), sseq.at(j).end(), sseq.at(j).begin(), ::toupper);
            alignments.at(j).append(sseq.at(j));
        }
    }

    if (mend.at(0) == sequences.at(0).size()-1) {
        for (int i = 0; i < reads.size(); i++) {
            if (reads.at(i).read_index(mend.at(0)) != sequences.at(i+1).size()-1) {
                extractSubsequences(mend.at(0), -1, sseq, false);
                // --- CREATE SUBSEQUENCE FASTA FILES ---
                writeFasta(directory + "temp/temp_sseq.fa", seq_names, sseq);

                // --- CALL MSA ALGORITHM ---
                system((directory + "mafft/mafft.bat " + directory + "temp/temp_sseq.fa > " + directory + "temp/temp_msa.fa").c_str());

                // --- READ MSA OUTPUT FILE ---
                sseq.clear();
                readFasta(directory + "temp/temp_msa.fa", sseq);
                for(int j = 0; j < sseq.size(); j++) {
                    transform(sseq.at(j).begin(), sseq.at(j).end(), sseq.at(j).begin(), ::toupper);
                    alignments.at(j).append(sseq.at(j));
                }
                break;
            }
        }
    } else {
        extractSubsequences(mend.at(0), -1, sseq, false);
        // --- CREATE SUBSEQUENCE FASTA FILES ---
        writeFasta(directory + "temp/temp_sseq.fa", seq_names, sseq);

        // --- CALL MSA ALGORITHM ---
        system((directory + "mafft/mafft.bat " + directory + "temp/temp_sseq.fa > " + directory + "temp/temp_msa.fa").c_str());

        // --- READ MSA OUTPUT FILE ---
        sseq.clear();
        readFasta(directory + "temp/temp_msa.fa", sseq);
        for(int j = 0; j < sseq.size(); j++) {
            transform(sseq.at(j).begin(), sseq.at(j).end(), sseq.at(j).begin(), ::toupper);
            alignments.at(j).append(sseq.at(j));
        }
    }

    // ---------- WRITE COMBINED ALIGNMENT RESULTS INTO OUTPUT FILE ----------

    writeFasta(res_file_name, seq_names, alignments);
}

int main(int argc, char* argv[])
{
    directory = argv[0];
    directory = directory.substr(0,directory.length()-5);
    ref_file_name = directory + "reference.fa";
    read_file_name = directory + "input.fa";
    res_file_name = directory + "alignment.fa";

    string argument;
    if (argc > 1) {
        for (int i = 1; i < argc; i++) {
            argument = argv[i];
            if (argument.at(0) == '-') {
                switch (argument.at(1))
                {
                case 's':
                    min_seed_size = stoi(argv[++i]);
                    break;
                case 'g':
                    max_merge_gap = stoi(argv[++i]);
                    break;
                case 'r':
                    ref_file_name = argv[++i];
                    break;
                case 'i':
                    read_file_name = argv[++i];
                    break;
                case 'o':
                    res_file_name = argv[++i];
                    break;
                default:
                    cerr << "[ERROR] INCORRECT ARGUMENTS" << endl;
                    return 1;
                }
            } else {
                cerr << "[ERROR] INCORRECT ARGUMENTS" << endl;
                return 1;
            }
        }
    }

    try {
        performMEMSA();
    } catch (int e) {
        cout << "calling MAFFT on unchanged input" << endl;
        performMafft();
    }

    return 0;
}
