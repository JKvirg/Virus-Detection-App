// CMSC 341 - Spring 26 - Project 4
#include "vdetect.h"
VDetect::VDetect(int size, hash_fn hash, prob_t probing = DEFPOLCY) {
    // Input Validation for Size
    if (size < MINPRIME) 
        m_currentCap = MINPRIME;
    else if (size > MAXPRIME) 
        m_currentCap = MAXPRIME;
    else {
        if (isPrime(size) == true) {
            m_currentCap = size;
        }
        else {
            m_currentCap = findNextPrime(size);
        }
    }

    // Initialize Member Variables
    m_hash = hash;
    m_currProbing = probing;

    m_currentSize = 0;
    m_currNumDeleted = 0;

    m_currentTable = new Virus* [m_currentCap];
    for (int i = 0; i < m_currentCap; i++) {
        m_currentTable[i] = nullptr;
    }

    // Initialize Member Variable for Old Table
    m_oldTable = nullptr;
    m_oldCap = 0;
    m_oldSize = 0;
    m_oldNumDeleted = 0;
    m_oldProbing = DEFPOLCY;

    m_transferIndex = 0;
    m_newPolicy = DEFPOLCY;
}

VDetect::~VDetect() {
    // Delete Current Table
    for (int i = 0; i < m_currentCap; i++) {
        delete m_currentTable[i];
    }
    delete[] m_currentTable;
    m_currentTable = nullptr;

    // Delete Old Table
    if (m_oldTable != nullptr) {
        for (int i = 0; i < m_oldCap; i++) {
            delete m_oldTable[i];
        }
        delete[] m_oldTable;
        m_oldTable = nullptr;
    }
}

void VDetect::changeProbPolicy(prob_t policy) {
    m_newPolicy = policy;
}

bool VDetect::insert(Virus virus) {
    // ID validation
    if (virus.getID() < MINID || virus.getID() > MAXID) {
        return false;
    }

    // Duplication Check
    if (getVirus(virus.getKey(), virus.getID()).getKey() != "") {
        return false;
    }

    // Transfer 25% if applicable
    if (m_oldTable != nullptr) {
        transferChunk();
    }

    // Insert into Hash Table
    if (insertHelper(m_currentTable, m_currentCap, m_currProbing, virus, m_currNumDeleted) == false) {
        return false;
    }
    m_currentSize++;

    // Rehash Conditions
    // Check Load Ratio
    // Check Deleted Ratio
    if (lambda() > 0.5 || deletedRatio() > 0.8) {
        if (m_oldTable == nullptr) {
            rehash();
        }
    }

    return true;
}

bool VDetect::remove(Virus virus) {
    int hashVal = m_hash(virus.getKey());
    int currIndex = hashVal % m_currentCap;
    int oldIndex = 0;
    if (m_oldCap != 0) {
        oldIndex = hashVal % m_oldCap;
    }
    int probe = 0;

    bool removed = false;
    // Search Current Table
    for (int i = 0; i < m_currentCap; i++) {
        // Get Index Value Based on Probing Policy
        // Case: Linear
        if (m_currProbing == LINEAR) {
            probe = (currIndex + i) % m_currentCap;
        }

        // Case: DoubleHash
        else if (m_currProbing == DOUBLEHASH) {
            probe = (currIndex + i * (11 - (hashVal % 11))) % m_currentCap;
        }

        // Case: Quadratic
        else {
            probe = (currIndex + (i * i)) % m_currentCap;
        }

        // Key Location is Empty
        if (m_currentTable[probe] == nullptr) break;

        // Key Location is Wanted Virus
        if (m_currentTable[probe]->getUsed() == true &&
            m_currentTable[probe]->getKey() == virus.getKey() &&
            m_currentTable[probe]->getID() == virus.getID()) {

            m_currentTable[probe]->setUsed(false);
            m_currNumDeleted++;
            removed = true;
            break;
        }
    }

    probe = 0;
    // Search Old Table
    if (m_oldTable != nullptr) {
        for (int i = 0; i < m_oldCap; i++) {
            // Get Index Value Based on Probing Policy
            // Case: Linear
            if (m_oldProbing == LINEAR) {
                probe = (oldIndex + i) % m_oldCap;
            }

            // Case: DoubleHash
            else if (m_oldProbing == DOUBLEHASH) {
                probe = (oldIndex + i * (11 - (hashVal % 11))) % m_oldCap;
            }

            // Case: Quadratic
            else {
                probe = (oldIndex + (i * i)) % m_oldCap;
            }

            if (m_oldTable[probe] == nullptr) break;

            // Key Location contains Wanted Virus
            if (m_oldTable[probe]->getUsed() == true &&
                m_oldTable[probe]->getKey() == virus.getKey() &&
                m_oldTable[probe]->getID() == virus.getID()) {

                m_oldTable[probe]->setUsed(false);
                m_oldNumDeleted++;
                removed = true;
                break;
            }
        }
    }

    // Transfer 25% Chunk if Applicable
    if (m_oldTable != nullptr) {
        transferChunk();
    }

    // Trigger Rehash if Applicable
    if (lambda() > 0.5 || deletedRatio() > 0.8) {
        if (m_oldTable == nullptr) {
            rehash();
        }
    }

    return removed;
}

const Virus VDetect::getVirus(string sequence, int id) const {
    // Sequence = key
    int hashVal = m_hash(sequence);
    int currIndex = hashVal % m_currentCap;
    int oldIndex = 0;
    if (m_oldCap != 0) {
        oldIndex = hashVal % m_oldCap;
    }
    int probe = 0;

    // Search Current Table
    for (int i = 0; i < m_currentCap; i++) {
        // Get Index Value Based on Probing Policy
        // Case: Linear
        if (m_currProbing == LINEAR) {
            probe = (currIndex + i) % m_currentCap;
        }

        // Case: DoubleHash
        else if (m_currProbing == DOUBLEHASH) {
            probe = (currIndex + i * (11 - (hashVal % 11))) % m_currentCap;
        }

        // Case: Quadratic
        else {
            probe = (currIndex + (i * i)) % m_currentCap;
        }

        // Key Location is Empty
        if (m_currentTable[probe] == nullptr) break;

        // Key Location contains Wanted Virus
        if (m_currentTable[probe]->getUsed() == true &&
            m_currentTable[probe]->getKey() == sequence &&
            m_currentTable[probe]->getID() == id) {
            
            return *m_currentTable[probe];
        }

    }

    probe = 0;
    // Search Old Table
    if (m_oldTable != nullptr) {
        for (int i = 0; i < m_oldCap; i++) {
            // Get Index Value Based on Probing Policy
            // Case: Linear
            if (m_oldProbing == LINEAR) {
                probe = (oldIndex + i) % m_oldCap;
            }

            // Case: DoubleHash
            else if (m_oldProbing == DOUBLEHASH) {
                probe = (oldIndex + i * (11 - (hashVal % 11))) % m_oldCap;
            }

            // Case: Quadratic
            else {
                probe = (oldIndex + (i * i)) % m_oldCap;
            }

            if (m_oldTable[probe] == nullptr) break;

            // Key Location contains Wanted Virus
            if (m_oldTable[probe]->getUsed() == true &&
                m_oldTable[probe]->getKey() == sequence &&
                m_oldTable[probe]->getID() == id) {

                return *m_oldTable[probe];
            }

        }
    }

    // Virus Doesn't Exist, Return Empty Virus Object
    return (Virus());
}

bool VDetect::updateID(Virus virus, int id) {
    Virus found = getVirus(virus.getKey(), virus.getID());

    // Virus Doesn't Exist
    if (found.getKey() == "") return false;

    // Remove Old Virus, Update ID, Insert Updateded Virus
    remove(virus);
    Virus updated = virus;
    updated.setID(id);
    insert(updated);

    return true;
}

float VDetect::lambda() const {
    return (static_cast<float>(m_currentSize) / m_currentCap);
}

float VDetect::deletedRatio() const {
    if (m_currentSize == 0) return 0;
    return (static_cast<float>(m_currNumDeleted) / m_currentSize);
}

void VDetect::dump() const {
    cout << "Dump for the current table: " << endl;
    if (m_currentTable != nullptr)
        for (int i = 0; i < m_currentCap; i++) {
            cout << "[" << i << "] : " << m_currentTable[i] << endl;
        }
    cout << "Dump for the old table: " << endl;
    if (m_oldTable != nullptr)
        for (int i = 0; i < m_oldCap; i++) {
            cout << "[" << i << "] : " << m_oldTable[i] << endl;
        }
}

bool VDetect::isPrime(int number) {
    bool result = true;
    for (int i = 2; i <= number / 2; ++i) {
        if (number % i == 0) {
            result = false;
            break;
        }
    }
    return result;
}

int VDetect::findNextPrime(int current) {
    //we always stay within the range [MINPRIME-MAXPRIME]
    //the smallest prime starts at MINPRIME
    if (current < MINPRIME) current = MINPRIME - 1;
    for (int i = current; i < MAXPRIME; i++) {
        for (int j = 2; j * j <= i; j++) {
            if (i % j == 0)
                break;
            else if (j + 1 > sqrt(i) && i != current) {
                return i;
            }
        }
    }
    //if a user tries to go over MAXPRIME
    return MAXPRIME;
}