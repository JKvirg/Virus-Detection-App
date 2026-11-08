// CMSC 341 - Spring 2026 - Project 4
#include "vdetect.h"
#include <random>
#include <vector>
enum RANDOM { UNIFORMINT, UNIFORMREAL, NORMAL };
class Random {
public:
    Random(int min, int max, RANDOM type = UNIFORMINT, int mean = 50, int stdev = 20) : m_min(min), m_max(max), m_type(type)
    {
        if (type == NORMAL) {
            //the case of NORMAL to generate integer numbers with normal distribution
            m_generator = std::mt19937(m_device());
            //the data set will have the mean of 50 (default) and standard deviation of 20 (default)
            //the mean and standard deviation can change by passing new values to constructor 
            m_normdist = std::normal_distribution<>(mean, stdev);
        }
        else if (type == UNIFORMINT) {
            //the case of UNIFORMINT to generate integer numbers
            // Using a fixed seed value generates always the same sequence
            // of pseudorandom numbers, e.g. reproducing scientific experiments
            // here it helps us with testing since the same sequence repeats
            m_generator = std::mt19937(10);// 10 is the fixed seed value
            m_unidist = std::uniform_int_distribution<>(min, max);
        }
        else { //the case of UNIFORMREAL to generate real numbers
            m_generator = std::mt19937(10);// 10 is the fixed seed value
            m_uniReal = std::uniform_real_distribution<double>((double)min, (double)max);
        }
    }
    void setSeed(int seedNum) {
        // we have set a default value for seed in constructor
        // we can change the seed by calling this function after constructor call
        // this gives us more randomness
        m_generator = std::mt19937(seedNum);
    }

    int getRandNum() {
        // this function returns integer numbers
        // the object must have been initialized to generate integers
        int result = 0;
        if (m_type == NORMAL) {
            //returns a random number in a set with normal distribution
            //we limit random numbers by the min and max values
            result = m_min - 1;
            while (result < m_min || result > m_max)
                result = m_normdist(m_generator);
        }
        else if (m_type == UNIFORMINT) {
            //this will generate a random number between min and max values
            result = m_unidist(m_generator);
        }
        return result;
    }

    double getRealRandNum() {
        // this function returns real numbers
        // the object must have been initialized to generate real numbers
        double result = m_uniReal(m_generator);
        // a trick to return numbers only with two deciaml points
        // for example if result is 15.0378, function returns 15.03
        // to round up we can use ceil function instead of floor
        result = std::floor(result * 100.0) / 100.0;
        return result;
    }

private:
    int m_min;
    int m_max;
    RANDOM m_type;
    std::random_device m_device;
    std::mt19937 m_generator;
    std::normal_distribution<> m_normdist;//normal distribution
    std::uniform_int_distribution<> m_unidist;//integer uniform distribution
    std::uniform_real_distribution<double> m_uniReal;//real uniform distribution

};

unsigned int hashCode(const string str);
string sequencer(int size, int seedNum);

// Testing/Debugging HashCode
unsigned int testHash(string str);

class Tester {
public:
    // Test Constructor
    bool testConstructorNonPrime() {
        VDetect detector(MINPRIME + 100, hashCode, QUADRATIC);

        if (detector.m_currentCap != 211) {
            cout << "FAILED: Capacity of Table is not the Expected Value.\n";
            return false;
        }
        if (detector.m_currProbing != QUADRATIC) {
            cout << "FAILED: Probing Policy is not the Expected Policy.\n";
            return false;
        }

        return true;
    }
    bool testConstructorPrime() {
        VDetect detector(MINPRIME + 110, hashCode, LINEAR);
        
        if (detector.m_currentCap != 211) {
            cout << "FAILED: Capacity of Table is not the Expected Value.\n";
            return false;
        }
        if (detector.m_currProbing != LINEAR) {
            cout << "FAILED: Probing Policy is not the Expected Policy.\n";
            return false;
        }

        return true;
    }
    bool testConstructorMaxError() {
        VDetect detector(MINPRIME - 1, hashCode, DOUBLEHASH);

        if (detector.m_currentCap != MINPRIME) {
            cout << "FAILED: Capacity of Table is not the Expected Value.\n";
            return false;
        }
        if (detector.m_currProbing != DOUBLEHASH) {
            cout << "FAILED: Probing Policy is not the Expected Policy.\n";
            return false;
        }

        return true;
    }
    bool testConstructorMinError() {
        VDetect detector(MAXPRIME + 1, hashCode, QUADRATIC);

        if (detector.m_currentCap != MAXPRIME) {
            cout << "FAILED: Capacity of Table is not the Expected Value.\n";
            return false;
        }
        if (detector.m_currProbing != QUADRATIC) {
            cout << "FAILED: Probing Policy is not the Expected Policy.\n";
            return false;
        }

        return true;
    }

    // Test Change Prob Policy
    bool testChangeProbPolicy() {
        VDetect detector(MINPRIME + 100, hashCode, QUADRATIC);

        detector.changeProbPolicy(LINEAR);

        if (detector.m_newPolicy != LINEAR || detector.m_currProbing != QUADRATIC) {
            cout << "FAILED: Current Probing has been changed or New Policy is not the Expected Policy.\n";
            return false;
        }

        return true;
    }

    // Test Lambda Function
    bool testLambdaFunction() {
        // Data Generation
        vector<Virus> dataList;
        Random RndID(MINID, MAXID);
        VDetect vdetect(MINPRIME, hashCode, DOUBLEHASH);
        
        for (int i = 0;i < 49;i++) {
            // generating random data
            Virus dataObj = Virus(sequencer(5, i), RndID.getRandNum(), true);
            // saving data for later use
            dataList.push_back(dataObj);
            // inserting data in to the VDetect object
            if (!vdetect.insert(dataObj)) cout << "Did not insert " << &dataObj << endl;
        }

        // Obtain Lambda Ratio Before Removing Nodes
        float ratio = vdetect.lambda();
        
        vdetect.remove(dataList[5]);
        vdetect.remove(dataList[15]);

        // Obtain Lambda Ratio After Removing Nodes
        float ratio2 = vdetect.lambda();

        if (ratio != ratio2) {
            cout << "FAILED: Lambda Ratio is not the expected value.\n";
            return false;
        }

        return true;
    }

    // Test Deleted Ratio Function
    bool testDeletedRatioFunction() {
        // Data Generation
        vector<Virus> dataList;
        Random RndID(MINID, MAXID);
        VDetect vdetect(MINPRIME, hashCode, DOUBLEHASH);

        for (int i = 0;i < 49;i++) {
            // generating random data
            Virus dataObj = Virus(sequencer(5, i), RndID.getRandNum(), true);
            // saving data for later use
            dataList.push_back(dataObj);
            // inserting data in to the VDetect object
            if (!vdetect.insert(dataObj)) cout << "Did not insert " << &dataObj << endl;
        }

        float ratio1 = vdetect.deletedRatio();

        vdetect.remove(dataList[5]);
        vdetect.remove(dataList[15]);
        vdetect.remove(dataList[23]);
        vdetect.remove(dataList[5]);
        vdetect.remove(dataList[2]);
        vdetect.remove(dataList[47]);

        float ratio2 = vdetect.deletedRatio();

        if (ratio1 != 0) {
            cout << "FAILED: Deleted Ratio before any Removals is not 0.\n";
            return false;
        }
        if (ratio2 == ratio1) {
            cout << "FAILED: Deleted Ratio after removals is same as before removals.\n";
            return false;
        }
        return true;
    }

    // Test Insert Function
    bool testInsert40() {
        vector<Virus> dataList;
        Random RndID(MINID, MAXID);
        VDetect vdetect(251, testHash, LINEAR);

        // Insert 10 Non Colliding Nodes
        for (int i = 0;i < 10;i++) {
            // generating random data
            Virus dataObj = Virus(to_string(i), MINID + i, true);
            // saving data for later use
            dataList.push_back(dataObj);
            // inserting data in to the VDetect object
            if (!vdetect.insert(dataObj)) cout << "Did not insert " << &dataObj << endl;
        }

        // Check Those 10 Non Colliding Nodes
        for (int i = 0;i < 10;i++) {
            if (*vdetect.m_currentTable[i] == dataList[i]) {
            }
            else {
                cout << "FAILED: Insertion of 10 non colliding nodes.\n";
                cout << "\tPosition: " << i << " is not the expected node.\n";
                return false;
            }
        }

        // Insert 10 Colliding Nodes with Linear Probing Type
        int j = 11;
        for (int i = 11;i < 21;i++) {
            // generating random data
            Virus dataObj = Virus(to_string(j), MINID + i, true);
            // saving data for later use
            dataList.push_back(dataObj);
            // inserting data in to the VDetect object
            if (!vdetect.insert(dataObj)) cout << "Did not insert " << &dataObj << endl;
        }

        // Check 10 Colliding Nodes with Linear Probing Type
        for (int i = 11;i < 21;i++) {
            if (*vdetect.m_currentTable[i] == dataList[i-1]) {
            }
            else {
                cout << "FAILED: Insertion of 10 colliding nodes of Linear.\n";
                cout << "\tPosition: " << i << " is not the expected node.\n";
                return false;
            }
        }

        // Insert 10 Colliding Nodes with DoubleHash Probing Type
        vdetect.m_currProbing = DOUBLEHASH;
        j = 25;
        for (int i = 25;i < 35 ;i++) {
            // generating random data
            Virus dataObj = Virus(to_string(j), MINID + i, true);
            // saving data for later use
            dataList.push_back(dataObj);
            // inserting data in to the VDetect object
            if (!vdetect.insert(dataObj)) cout << "Did not insert " << &dataObj << endl;
        }
        // Check 10 Colliding Nodes with DoubleHash Probing Type
        for (int i = 0;i < 10;i++) {
            if (*vdetect.m_currentTable[25+ (8*i)] == dataList[i + 20]) {
            }
            else {
                cout << "FAILED: Insertion of 10 colliding nodes of Linear.\n";
                cout << "\tPosition: " << (25 + (8 * i)) << " is not the expected node.\n";
                return false;
            }
        }

        // Insert 10 Colliding Nodes with Quadratic Probing Type
        vdetect.m_currProbing = QUADRATIC;
        j = 100;
        for (int i = 100;i < 110;i++) {
            // generating random data
            Virus dataObj = Virus(to_string(j), MINID + i, true);
            // saving data for later use
            dataList.push_back(dataObj);
            // inserting data in to the VDetect object
            if (!vdetect.insert(dataObj)) cout << "Did not insert " << &dataObj << endl;
        }
        // Check 10 Colliding Nodes with Quadratic Probing Type
        for (int i = 0;i < 10;i++) {
            if (*vdetect.m_currentTable[((100 +(i*i)) % 251)] == dataList[i + 30]) {
            }
            else {
                cout << "FAILED: Insertion of 10 colliding nodes of Linear.\n";
                cout << "\tPosition: " << ((100 + (i * i)) % 251) << " is not the expected node.\n";
                return false;
            }
        }

        if (vdetect.m_currentSize != 40) {
            cout << "FAILED: Size Tracker is not equal to 40 after inserting 40 nodes.\n";
            return false;
        }

        return true;
    }
    
    // Test Remove Function
    bool testRemove20() {
        // Data Generation
        vector<Virus> dataList;
        Random RndID(MINID, MAXID);
        VDetect vdetect(MINPRIME, hashCode, QUADRATIC);

        for (int i = 0;i < 49;i++) {
            // generating random data
            Virus dataObj = Virus(sequencer(5, i), RndID.getRandNum(), true);
            // saving data for later use
            dataList.push_back(dataObj);
            // inserting data in to the VDetect object
            if (!vdetect.insert(dataObj)) cout << "Did not insert " << &dataObj << endl;
        }

        for (int i = 0; i < 20; i++) {
            if (vdetect.remove(dataList[i]) == false) {
                cout << "FAILED: Did not remove requested node inside the Hash Table.\n";
                return false;
            }
        }

        if (vdetect.m_currentSize != 49) {
            cout << "FAILED: Size of the Hash Table is not the expected Value.\n";
            return false;
        }

        if (vdetect.m_currNumDeleted != 20) {
            cout << "FAILED: Current Num of Deleted is not the expected Value.\n";
            return false;
        }

        return true;
    }
    
    // Test GetVirus Function
    bool testGetVirusNormal() {
        // Data Generation
        vector<Virus> dataList;
        Random RndID(MINID, MAXID);
        VDetect vdetect(MINPRIME, hashCode, QUADRATIC);

        for (int i = 0;i < 49;i++) {
            // generating random data
            Virus dataObj = Virus(sequencer(5, i), RndID.getRandNum(), true);
            // saving data for later use
            dataList.push_back(dataObj);
            // inserting data in to the VDetect object
            if (!vdetect.insert(dataObj)) cout << "Did not insert " << &dataObj << endl;
        }

        Virus validation = vdetect.getVirus(dataList[10].getKey(), dataList[10].getID());

        if (validation.getKey() != dataList[10].getKey() || validation.getID() != dataList[10].getID()) {
            cout << "FAILED: Returned Virus is not the Same.\n";
            return false;
        }
        return true;
    }
    bool testGetVirusError() {
        // Data Generation
        vector<Virus> dataList;
        Random RndID(MINID + 1, MAXID);
        VDetect vdetect(MINPRIME, hashCode, QUADRATIC);

        for (int i = 0;i < 49;i++) {
            // generating random data
            Virus dataObj = Virus(sequencer(5, i), RndID.getRandNum(), true);
            // saving data for later use
            dataList.push_back(dataObj);
            // inserting data in to the VDetect object
            if (!vdetect.insert(dataObj)) cout << "Did not insert " << &dataObj << endl;
        }
        Virus error("AAAAA", MINID, true);
        Virus validation = vdetect.getVirus(error.getKey(), error.getID());

        if (validation.getKey() != "" || validation.getID() != 0) {
            cout << "FAILED: Returned Virus is not an Empty Virus Object.\n";
            return false;
        }
        return true;
    }
    bool testGetVirusNonCollide() {
        // Data Generation
        vector<Virus> dataList;
        Random RndID(MINID, MAXID);
        VDetect vdetect(MINPRIME, hashCode, QUADRATIC);

        for (int i = 0;i < 49;i++) {
            // generating random data
            Virus dataObj = Virus(sequencer(5, i), RndID.getRandNum(), true);
            // saving data for later use
            dataList.push_back(dataObj);
            // inserting data in to the VDetect object
            if (!vdetect.insert(dataObj)) cout << "Did not insert " << &dataObj << endl;
        }

        
        for (int i = 48; i >= 0; i--) {
            Virus validation = vdetect.getVirus(dataList[i].getKey(), dataList[i].getID());

            if (validation.getKey() != dataList[i].getKey() || validation.getID() != dataList[i].getID()) {
                cout << "FAILED: Returned Object is not the expected Object.\n";
                return false;
            }


        }

        return true;
    }
    bool testGetVirusCollide() {
        // Data Generation
        vector<Virus> dataList;
        VDetect vdetect(MINPRIME, hashCode, QUADRATIC);

        for (int i = 0;i < 49;i++) {
            // generating random data
            Virus dataObj = Virus("AAAAA", MINID + i, true);
            // saving data for later use
            dataList.push_back(dataObj);
            // inserting data in to the VDetect object
            if (!vdetect.insert(dataObj)) cout << "Did not insert " << &dataObj << endl;
        }


        for (int i = 48; i >= 0; i--) {
            Virus validation = vdetect.getVirus(dataList[i].getKey(), dataList[i].getID());

            if (validation.getKey() != dataList[i].getKey() || validation.getID() != dataList[i].getID()) {
                cout << "FAILED: Returned Object is not the expected Object.\n";
                return false;
            }


        }

        return true;
    }
    
    // Test UpdateID Function
    bool testUpdateID() {
        // Data Generation
        vector<Virus> dataList;
        Random RndID(MINID, MAXID - 1);
        VDetect vdetect(MINPRIME, hashCode, QUADRATIC);

        for (int i = 0;i < 49;i++) {
            // generating random data
            Virus dataObj = Virus(sequencer(5, i), RndID.getRandNum(), true);
            // saving data for later use
            dataList.push_back(dataObj);
            // inserting data in to the VDetect object
            if (!vdetect.insert(dataObj)) cout << "Did not insert " << &dataObj << endl;
        }

        vdetect.updateID(dataList[10], MAXID);

        Virus validation = vdetect.getVirus(dataList[10].getKey(), MAXID);

        if (validation.getKey() != dataList[10].getKey() || validation.getID() != MAXID) {
            cout << "FAILED: Id of returned Virus is not the expected Value.\n";
            return false;
        }
        
        return true;
    }

private:
};

int main() {
    Tester tester;
    // Constructor Tests
    if (tester.testConstructorNonPrime() == true) cout << "PASSED: Test of Constructor NON PRIME.\n";
    if (tester.testConstructorPrime() == true) cout << "PASSED: Test of Constructor PRIME.\n";
    if (tester.testConstructorMaxError() == true) cout << "PASSED: Test of Constructor MAX Error.\n";
    if (tester.testConstructorMinError() == true) cout << "PASSED: Test of Constructor MIN Error.\n";
    // Change Prob Policy Test
    if (tester.testChangeProbPolicy() == true) cout << "PASSED: Test of Change Prob Policy.\n";
    // Lambda Func Test
    if (tester.testLambdaFunction() == true) cout << "PASSED: Test of Lambda Function.\n";
    // DeletedRatio Func Test
    if (tester.testDeletedRatioFunction() == true) cout << "PASSED: Test of Deleted Ratio Function.\n";
    // Insert Func Tests
    if (tester.testInsert40() == true) cout << "PASSED: Test of Insert 40.\n";
    // Remove Func Tests
    if (tester.testRemove20() == true) cout << "PASSED: Test of Remove 40.\n";
    // GetVirus Func Tests
    if (tester.testGetVirusNormal() == true) cout << "PASSED: Test of Get Virus Normal.\n";
    if (tester.testGetVirusError() == true) cout << "PASSED: Test of Get Virus Error.\n";
    if (tester.testGetVirusNonCollide() == true) cout << "PASSED: Test of Get Virus Non-Collide.\n";
    if (tester.testGetVirusCollide() == true) cout << "PASSED: Test of Get Virus Collide.\n";
    // UpdateID Func Tests
    if (tester.testUpdateID() == true) cout << "PASSED: Test of UpdateID.\n";

    return 0;
}
unsigned int hashCode(const string str) {
    unsigned int val = 0;
    const unsigned int thirtyThree = 33;
    for (int i = 0; i < str.length(); i++)
        val = val * thirtyThree + str[i];
    return val;
}
string sequencer(int size, int seedNum) {
    //this function returns a random DNA sequence
    string sequence = "";
    Random rndObject(0, 3);
    rndObject.setSeed(seedNum);
    for (int i = 0;i < size;i++) {
        sequence = sequence + ALPHA[rndObject.getRandNum()];
    }
    return sequence;
}

// Used for Testing/Debugging
// Takes str in form of "123"
// Returns as an integer of 123
unsigned int testHash(string str) {
    unsigned int val = 0;
    val = stoi(str);
    return val;
}
