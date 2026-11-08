// CMSC 341 - Spring 2026 - Project 4
#ifndef VDETECT_H
#define VDETECT_H
#include <iostream>
#include <string>
#include "math.h"
using namespace std;
class Grader;   // forward declaration, will be used for grdaing
class Tester;   // forward declaration, will be used for testing
class Virus;    // forward declaration
class VDetect;  // forward declaration
const int MINID = 100000;
const int MAXID = 999999;
const int MINPRIME = 101;   // Min size for hash table
const int MAXPRIME = 99991; // Max size for hash table
#define EMPTY Virus("",0)
#define DELETED Virus("DELETED")
#define DELETEDKEY "DELETED"
typedef unsigned int (*hash_fn)(string); // declaration of hash function
enum prob_t { QUADRATIC, DOUBLEHASH, LINEAR }; // types of collision handling policy
#define DEFPOLCY QUADRATIC
// from https://www.iarcs.org.in/inoi/online-study-material/problems/dna-soln.php
//  Biologically, a strand of DNA is a chain of nucleotides 
//  Adenine, Cytosine, Guanine, and Thymine. The four nucleotides 
//  are represented by characters A, C, G, and T, respectively. 
//  Thus, a strand of DNA can be represented by a string of 
//  these four characters. We call such a string a DNA sequence.
const int MAX = 4;
const char ALPHA[MAX] = { 'A', 'C', 'G', 'T' };

class Virus {
public:
    friend class Grader;
    friend class Tester;
    friend class VDetect;
    Virus(string key = "", int id = 0, bool used = false) {
        m_key = key; m_id = id; m_used = used;
    }
    string getKey() const { return m_key; }
    int getID() const { return m_id; }
    void setKey(string key) { m_key = key; }
    void setID(int id) { m_id = id; }
    bool getUsed() const { return m_used; }
    void setUsed(bool used) { m_used = used; }
    const Virus& operator=(const Virus& rhs) {
        if (this != &rhs) {
            m_key = rhs.m_key;
            m_id = rhs.m_id;
            m_used = rhs.m_used;
        }
        return *this;
    }
    // the following function is a friend function
    friend ostream& operator<<(ostream& sout, const Virus* virus) {
        if ((virus != nullptr) && !(virus->getKey().empty()))
            sout << virus->getKey() << " (" << virus->getID() << ", " << virus->getUsed() << ")";
        else
            sout << "";
        return sout;
    }
    // the following function is a friend function
    friend bool operator==(const Virus& lhs, const Virus& rhs) {
        // since the uniqueness of an object is defined by key and ID
        // the equality operator considers only those two criteria
        return ((lhs.getKey() == rhs.getKey()) && (lhs.getID() == rhs.getID()));
    }
    // the following function is a class function
    bool operator==(const Virus*& rhs) {
        // since the uniqueness of an object is defined by sequence and ID
        // the equality operator considers only those two criteria
        return ((getKey() == rhs->getKey()) && (getID() == rhs->getID()));
    }
private:
    string m_key;   // the search string used as key in the hash table
    int m_id;       // a unique ID number identifying the object
    // the following variable is used for lazy delete scheme in hash table
    // if it is set to false, it means the bucket in the hash table is free for insert
    // if it is set to true, it means the bucket contains live data, and we cannot overwrite it
    bool m_used;
};

class VDetect {
public:
    friend class Grader;
    friend class Tester;
    VDetect(int size, hash_fn hash, prob_t probing);
    ~VDetect();
    // Returns Load factor of the new table
    float lambda() const;
    // Returns the ratio of deleted slots in the new table
    float deletedRatio() const;
    // insert only happens in the new table
    bool insert(Virus virus);
    // remove can happen from either table
    bool remove(Virus virus);
    // find can happen in either table
    const Virus getVirus(string key, int id) const;
    // update the information
    bool updateID(Virus virus, int id);
    void changeProbPolicy(prob_t policy);
    void dump() const;

private:
    hash_fn    m_hash;          // hash function
    prob_t     m_newPolicy;     // stores the change of policy request

    Virus** m_currentTable;  // hash table
    int        m_currentCap;    // hash table size (capacity)
    int        m_currentSize;   // current number of entries
    // m_currentSize includes deleted entries 
    int        m_currNumDeleted;// number of deleted entries
    prob_t     m_currProbing;       // collision handling policy

    Virus** m_oldTable;      // hash table
    int        m_oldCap;        // hash table size (capacity)
    int        m_oldSize;       // current number of entries
    // m_oldSize includes deleted entries
    int        m_oldNumDeleted; // number of deleted entries
    prob_t     m_oldProbing;    // collision handling policy

    int        m_transferIndex; // this can be used as a temporary place holder
    // during incremental transfer to scanning the table

//private helper functions
    bool isPrime(int number);
    int findNextPrime(int current);

    // Inserts Virus into Current Table
    bool insertHelper(Virus** table, int cap, prob_t policy, Virus virus, int& numOfDeleted) {
        int hashVal = m_hash(virus.getKey());
        int index = hashVal % cap;
        
        // Insert
        int probe = 0;
        int firstDeleted = -1; // Used for Insertion if no duplicate Exist
        for (int i = 0; i < cap; i++) {
            // Get Index Value Based on Probing Policy
            // Case: Linear
            if (policy == LINEAR) {
                probe = (index + i) % cap;
            }

            // Case: DoubleHash
            else if (policy == DOUBLEHASH) {
                probe = (index + i * (11 - (hashVal % 11))) % cap;
            }

            // Case: Quadratic
            else {
                probe = (index + (i * i)) % cap;
            }

            if (table[probe] == nullptr) {
                if (firstDeleted != -1) {
                    probe = firstDeleted;
                }

                table[probe] = new Virus(virus);
                table[probe]->setUsed(true);
                return true;
            }

            // Check Duplicate
            if (table[probe]->getUsed() == true &&
                table[probe]->getKey() == virus.getKey() &&
                table[probe]->getID() == virus.getID()) {
                
                return false;
            }

            if (table[probe]->getUsed() == false && firstDeleted == -1) {
                firstDeleted = probe;
            }
        }

        if (firstDeleted != -1) {
            table[firstDeleted] = new Virus(virus);
            table[firstDeleted]->setUsed(true);
            numOfDeleted--;
            return true;
        }
        return false;
    }

    // Starts Rehash Process
    void rehash() {
        // Move Current Table to Old Variables
        m_oldTable = m_currentTable;
        m_oldCap = m_currentCap;
        m_oldSize = m_currentSize;
        m_oldNumDeleted = m_currNumDeleted;
        m_oldProbing = m_currProbing;
        
        // Calculate New Capacity
        int newCap = findNextPrime((m_oldSize - m_oldNumDeleted) * 4);

        // Assign New Prob Policy
        m_currProbing = m_newPolicy;

        // Create New Hash Table
        m_currentCap = newCap;
        m_currentTable = new Virus * [m_currentCap];
        for (int i = 0; i < m_currentCap; i++) {
            m_currentTable[i] = nullptr;
        }

        m_currentSize = 0;
        m_currNumDeleted = 0;
        m_transferIndex = 0;
    }

    // Transfer 25% of the Old Hash Table to Curr Hash Table
    void transferChunk() {
        if (m_oldTable == nullptr) return;

        // Obtain Integer Value containing 25% of the Table
        int chunk = m_oldCap / 4;

        // Transfer Virus Nodes
        for (int i = 0; i < chunk && m_transferIndex < m_oldCap; i++) {
            Virus* v = m_oldTable[m_transferIndex];

            // Check Position isn't Nullptr
            if (v != nullptr) {
                // Check if Virus is Live
                if (v->getUsed() == true) {
                    if (insertHelper(m_currentTable, m_currentCap, m_currProbing, *v, m_currNumDeleted) == true) {
                        m_currentSize++;
                        m_oldSize--;
                    }
                }
                else {
                    m_oldSize--;
                    m_oldNumDeleted--;
                }

                delete v;
                m_oldTable[m_transferIndex] = nullptr;
            }

            m_transferIndex++;
        }

        // Done Transfer
        if (m_transferIndex >= m_oldCap) {
            delete[] m_oldTable;
            m_oldTable = nullptr;
        }
    }
};
#endif