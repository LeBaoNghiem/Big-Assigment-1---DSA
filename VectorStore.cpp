#include "VectorStore.h"

// ----------------- ArrayList Implementation -----------------

template <class T>
ArrayList<T>::ArrayList(int initCapacity) : capacity(initCapacity), count(0)
{
    if (capacity <= 0)
        data = nullptr;
    else
        data = new T[initCapacity];
}

template <class T>
ArrayList<T>::ArrayList(const ArrayList<T> &other)
{
    capacity = other.capacity;
    count = other.count;
    data = new T[capacity];

    for (int i = 0; i < count; ++i)
        data[i] = other.data[i];
}

template <class T>
ArrayList<T>::~ArrayList()
{
    delete[] data;
}

template <class T>
ArrayList<T> &ArrayList<T>::operator=(const ArrayList<T> &other)
{
    if (this == &other)
        return *this;

    delete[] data;

    capacity = other.capacity;
    count = other.count;
    data = new T[capacity];

    for (int i = 0; i < count; ++i)
        data[i] = other.data[i];

    return *this;
}

template <class T>
void ArrayList<T>::ensureCapacity(int cap)
{
    if (cap > capacity)
    {
        capacity = (capacity) ? capacity * 3 / 2 : cap;
        T *newData = new T[capacity];

        for (int i = 0; i < count - 1; ++i)
        {
            newData[i] = data[i];
        }

        delete[] data;
        data = newData;
    }
}

template <class T>
void ArrayList<T>::add(T e)
{
    ensureCapacity(++count);
    data[count - 1] = e;
}

template <class T>
void ArrayList<T>::add(int index, T e)
{
    if (index < 0 || index > count)
        throw out_of_range("Index is invalid!");

    ensureCapacity(++count);
    for (int i = count - 1; i > index; --i)
        data[i] = data[i - 1];

    data[index] = e;
}

template <class T>
T ArrayList<T>::removeAt(int index)
{
    if (index < 0 || index >= count)
        throw out_of_range("Index is invalid!");

    T removedEle = data[index];
    for (int i = index; i < count - 1; ++i)
        data[i] = data[i + 1];

    --count;
    return removedEle;
}

template <class T>
bool ArrayList<T>::empty() const
{
    return count == 0;
}

template <class T>
int ArrayList<T>::size() const
{
    return count;
}

template <class T>
void ArrayList<T>::clear()
{
    delete[] data;
    T *newData = new T[10];
    data = newData;
    capacity = 10;
    count = 0;
}

template <class T>
T &ArrayList<T>::get(int index) const
{
    if (index < 0 || index >= count)
        throw out_of_range("Index is invalid!");

    return data[index];
}

template <class T>
void ArrayList<T>::set(int index, T e)
{
    if (index < 0 || index >= count)
        throw out_of_range("Index is invalid!");

    data[index] = e;
}

template <class T>
int ArrayList<T>::indexOf(T item) const
{
    for (int i = 0; i < count; ++i)
        if (data[i] == item)
            return i;

    return -1;
}

template <class T>
bool ArrayList<T>::contains(T item) const
{
    return indexOf(item) != -1;
}

template <class T>
string ArrayList<T>::toString(string (*item2str)(T &)) const
{
    ostringstream oss;

    oss << "[";
    for (int i = 0; i < count - 1; ++i)
    {
        if (item2str)
            oss << item2str(data[i]);
        else
            oss << data[i];
        oss << ", ";
    }
    if (count >= 1)
        if (item2str) oss << item2str(data[count - 1]);
        else oss <<  data[count - 1];
    oss << "]";

    return oss.str();
}

template <class T>
typename ArrayList<T>::Iterator ArrayList<T>::begin()
{
    return Iterator(this, 0);
}

template <class T>
typename ArrayList<T>::Iterator ArrayList<T>::end()
{
    return Iterator(this, count);
}

// ----------------- Iterator of ArrayList Implementation -----------------
template <class T>
ArrayList<T>::Iterator::Iterator(ArrayList<T> *pList, int index) : pList(pList)
{
    // TODO
    if (index < 0 || index > pList->count)
        throw out_of_range("Index is invalid!");

    cursor = index;
}

template <class T>
typename ArrayList<T>::Iterator &ArrayList<T>::Iterator::operator=(const Iterator &other)
{
    pList = other.pList;
    cursor = other.cursor;

    return *this;
}

template <class T>
T &ArrayList<T>::Iterator::operator*()
{
    if (cursor < 0 || cursor >= pList->count)
        throw out_of_range("Iterator is out of range!");

    return pList->data[cursor];
}

template <class T>
bool ArrayList<T>::Iterator::operator!=(const Iterator &other) const
{
    return cursor != other.cursor || pList != other.pList;
}

template <class T>
typename ArrayList<T>::Iterator &ArrayList<T>::Iterator::operator++()
{
    if (cursor >= pList->count)
        throw out_of_range("Iterator cannot advance past end!");

    ++cursor;
    return *this;
}

template <class T>
typename ArrayList<T>::Iterator ArrayList<T>::Iterator::operator++(int)
{
    if (cursor >= pList->count)
        throw out_of_range("Iterator cannot advance past end!");

    Iterator temp = *this;
    ++cursor;
    return temp;
}

template <class T>
typename ArrayList<T>::Iterator &ArrayList<T>::Iterator::operator--()
{
    if (cursor <= 0)
        throw out_of_range("Iterator cannot move before begin!");

    --cursor;
    return *this;
}

template <class T>
typename ArrayList<T>::Iterator ArrayList<T>::Iterator::operator--(int)
{
    if (cursor <= 0)
        throw out_of_range("Iterator cannot move before begin!");

    Iterator temp = *this;
    --cursor;
    return temp;
}

// ----------------- SinglyLinkedList Implementation -----------------
template <class T>
SinglyLinkedList<T>::SinglyLinkedList() : head(nullptr), tail(nullptr), count(0) {}

template <class T>
SinglyLinkedList<T>::~SinglyLinkedList() {
    clear();
}

template <class T>
void SinglyLinkedList<T>::add(T e) {
    Node* newNode = new Node(e);

    if (tail) tail->next = newNode;
    else {
        head = newNode;
    }

    tail = newNode;
    count++;
}

template <class T>
void SinglyLinkedList<T>::add(int index, T e) {
    if (index < 0 || index > count) throw out_of_range("Index is invalid!");

    int currentIndex = 0;
    Node* current = head;
    Node* previous = nullptr;
    Node* newNode = new Node(e);

    while (currentIndex < index) {
        previous = current;
        current = current->next;
        currentIndex++;
    }

    if (!previous) {
        newNode->next = head;
        head = newNode;
        if (!tail) tail = newNode;

    } else {
        previous->next = newNode;
        newNode->next = current;
        if (!current) tail = newNode;
    }
    
    count++;
}

template <class T>
T SinglyLinkedList<T>::removeAt(int index) {
    if (index < 0 || index >= count) throw out_of_range("Index is invalid!");

    int currentIndex = 0;
    Node* current = head;
    Node* previous = nullptr;

    while (currentIndex < index) {
        previous = current;
        current = current->next;
        currentIndex++;
    }

    T removedData = current->data;

    if (!previous) {
        head = head->next;
        if (!head) tail = nullptr;

    } else {
        previous->next = current->next;
        if (!current->next) tail = previous;
    }

    delete current;
    count--;
    return removedData;
}

template <class T>
bool SinglyLinkedList<T>::removeItem(T item) {
    Node* current = head;
    Node* previous = nullptr;

    while (current) {
        if (current->data == item) {
            if (!previous) {
                head = head->next;
                if (!head) tail = nullptr;

            } else {
                previous->next = current->next;
                if (!current->next) tail = previous;
            }

            delete current;
            count--;
            return true;
        }

        previous = current;
        current = current->next;
    }

    return false;
}

template <class T>
bool SinglyLinkedList<T>::empty() const {
    return count == 0;
}

template <class T>
int SinglyLinkedList<T>::size() const {
    return count;
}

template <class T>
void SinglyLinkedList<T>::clear() {
    Node* current = head;

    while (current) {
        Node* nextNode = current->next;
        delete current;
        current = nextNode;
    }

    head = nullptr;
    tail = nullptr;
    count = 0;
}

template <class T>
T& SinglyLinkedList<T>::get(int index) const {
    if (index < 0 || index >= count) throw out_of_range("Index is invalid!");

    int currentIndex = 0;
    Node* current = head;

    while (currentIndex < index) {
        current = current->next;
        currentIndex++;
    }

    return current->data;
}

template <class T>
int SinglyLinkedList<T>::indexOf(T item) const {
    int currentIndex = 0;
    Node* current = head;

    while (current) {
        if (current->data == item) return currentIndex;
        current = current->next;
        currentIndex++;
    }

    return -1;
}

template <class T>
bool SinglyLinkedList<T>::contains(T item) const {
    if (indexOf(item) != -1) return true;
    return false;
}

template <class T>
string SinglyLinkedList<T>::toString(string (*item2str)(T &)) const {
    ostringstream oss;
    Node* current = head;

    while (current) {
        oss << '[';
        if (item2str) oss << item2str(current->data);
        else oss << current->data;
        oss << ']';


        if (current->next) oss << "->";
        current = current->next;
    }

    return oss.str();
}

template <class T>
typename SinglyLinkedList<T>::Iterator SinglyLinkedList<T>::begin() {
    return Iterator(head);
}

template <class T>
typename SinglyLinkedList<T>::Iterator SinglyLinkedList<T>::end() {
    return Iterator(nullptr);
}


// ----------------- Iterator of SinglyLinkedList Implementation -----------------
template <class T>
SinglyLinkedList<T>::Iterator::Iterator(Node *node) : current(node) {}

template <class T>
typename SinglyLinkedList<T>::Iterator& SinglyLinkedList<T>::Iterator::operator=(const Iterator &other) {
    this->current = other.current;
    return *this;
}

template <class T>
T& SinglyLinkedList<T>::Iterator::operator*() {
    if (!current) throw out_of_range("Iterator is out of range!");

    return current->data;
}

template <class T>
bool SinglyLinkedList<T>::Iterator::operator!=(const Iterator &other) const {
    return this->current != other.current;
}

template <class T>
typename SinglyLinkedList<T>::Iterator& SinglyLinkedList<T>::Iterator::operator++() {
    if (!current) throw out_of_range("Iterator cannot advance past end!");
    
    current = current->next;
    return *this;
}

template <class T>
typename SinglyLinkedList<T>::Iterator SinglyLinkedList<T>::Iterator::operator++(int) {
    if (!current) throw out_of_range("Iterator cannot advance past end!");

    Iterator temp = *this;
    current = current->next;
    return temp;
}


// ----------------- VectorStore Implementation -----------------

VectorStore::VectorStore(int dimension, EmbedFn embeddingFunction) 
    : dimension(dimension), count(0), embeddingFunction(embeddingFunction) {}

VectorStore::~VectorStore() {
    clear();
}

VectorStore::VectorRecord::VectorRecord(int id, const string &rawText, SinglyLinkedList<float> *vector)
    : id(id), rawText(rawText), rawLength(rawText.length()), vector(vector) {}

VectorStore::VectorRecord::~VectorRecord() {
    delete vector;
}

int VectorStore::size() const {
    return count;
}

bool VectorStore::empty() const {
    return count == 0;
}

void VectorStore::clear() {
    for (int i = 0; i < count; ++i) delete records.get(i);
    records.clear();

    count = 0;
}

SinglyLinkedList<float>* VectorStore::preprocessing(string rawText) {
    if (embeddingFunction == nullptr) return nullptr;
    SinglyLinkedList<float>* vector = embeddingFunction(rawText);

    if (vector->size() < dimension) {
        int padding = dimension - vector->size();
        for (int i = 0; i < padding; ++i) vector->add(0.0f);

    } else if (vector->size() > dimension) while (vector->size() > dimension) vector->removeAt(vector->size() - 1);

    return vector;
}

void VectorStore::addText(string rawText) {
    SinglyLinkedList<float>* vector = preprocessing(rawText);

    int nextID;
    if (!records.empty()) nextID = records.get(count - 1)->id + 1;
    else nextID = 1;

    VectorRecord* newRecord = new VectorRecord(nextID, rawText, vector);
    records.add(newRecord);
    count++;
}

void VectorStore::checkIndex(int index) const {
    if (index < 0 || index >= count) throw out_of_range("Index is invalid!");
}

SinglyLinkedList<float>& VectorStore::getVector(int index) {
    checkIndex(index);

    return *(records.get(index)->vector);
}

string VectorStore::getRawText(int index) const {
    checkIndex(index);

    return records.get(index)->rawText;
}

int VectorStore::getId(int index) const {
    checkIndex(index);

    return records.get(index)->id;
}

bool VectorStore::removeAt(int index) {
    checkIndex(index);

    delete records.get(index);
    records.removeAt(index);

    count--;

    return true;
}

bool VectorStore::updateText(int index, string newRawText) {
    checkIndex(index);

    SinglyLinkedList<float>* newVector = preprocessing(newRawText);
    VectorRecord* record = records.get(index);

    delete record->vector;
    record->vector = newVector;

    record->rawText = newRawText;
    record->rawLength = newRawText.length();

    return true;
}

void VectorStore::setEmbeddingFunction(EmbedFn newEmbeddingFunction) {
    embeddingFunction = newEmbeddingFunction;
}

void VectorStore::forEach(void (*action)(SinglyLinkedList<float> &, int, string &)) {
    for (int i = 0; i < count; ++i) {
        VectorRecord* record = records.get(i);
        action(*(record->vector), record->rawLength, record->rawText);
    }
}

double VectorStore::cosineSimilarity(const SinglyLinkedList<float> &v1, const SinglyLinkedList<float> &v2) const {
    if (!v1.size() && !v2.size()) return 0.0;

    double dotProduct = 0.0;
    double lengthV1 = 0.0;
    double lengthV2 = 0.0;

    for (int i = 0; i < v1.size(); ++i) {
        double val1 = v1.get(i);
        double val2 = v2.get(i);

        dotProduct += val1 * val2;
        lengthV1 += val1 * val1;
        lengthV2 += val2 * val2;
    }

    if (fabs(lengthV1) < EPSILON || fabs(lengthV2) < EPSILON) return 0.0;

    return dotProduct / (sqrt(lengthV1) * sqrt(lengthV2));
}

double VectorStore::l1Distance(const SinglyLinkedList<float> &v1, const SinglyLinkedList<float> &v2) const {
    if (!v1.size() && !v2.size()) return 0.0;

    double distance = 0.0;

    for (int i = 0; i < v1.size(); ++i) 
        distance += fabs(v1.get(i) - v2.get(i));

    return distance;
}

double VectorStore::l2Distance(const SinglyLinkedList<float> &v1, const SinglyLinkedList<float> &v2) const {
    if (!v1.size() && !v2.size()) return 0.0;

    double distance = 0.0;

    for (int i = 0; i < v1.size(); ++i) {
        double diff = v1.get(i) - v2.get(i);
        distance += diff * diff;
    }

    return sqrt(distance);
}

int VectorStore::caculatingFirst(const SinglyLinkedList<float> &query, double (VectorStore::*func)(const SinglyLinkedList<float> &, const SinglyLinkedList<float> &) const, bool minimize) const {
    int bestIndex = 0;
    double bestValue = minimize ? 1e308 : -1e308;

    for (int i = 0; i < count; ++i) {
        double value = (this->*func)(query, *(records.get(i)->vector));
        if ((minimize && value + EPSILON < bestValue ) || (!minimize && value > bestValue + EPSILON)) {
            bestValue = value;
            bestIndex = i;
        }
    }

    return bestIndex;
}

void VectorStore::checkMetric(const string &metric) const {
    if (metric != "cosine" && metric != "euclidean" && metric != "manhattan") throw invalid_metric();
}

int VectorStore::findNearest(const SinglyLinkedList<float> &query, const string &metric) const {
    checkMetric(metric);

    if (metric == "cosine") return caculatingFirst(query, &VectorStore::cosineSimilarity, false);
    else if (metric == "euclidean") return caculatingFirst(query, &VectorStore::l2Distance, true);
    else return caculatingFirst(query, &VectorStore::l1Distance, true);
}

void VectorStore::merge(Data arr[], Data left[], int leftSize, Data right[], int rightSize, bool minimize) const {
    int i = 0, j = 0, k = 0;

    if (minimize) {
        while (j < leftSize && k < rightSize) {
            if (left[j].distance < right[k].distance - EPSILON) {
                arr[i++] = left[j++];

            } else if (right[k].distance < left[j].distance - EPSILON) {
                 arr[i++] = right[k++];

            } else {
                if (left[j].index < right[k].index) arr[i++] = left[j++];
                else arr[i++] = right[k++];
            }
        }

    } else {
        while (j < leftSize && k < rightSize) {
            if (left[j].distance > right[k].distance + EPSILON) {
                arr[i++] = left[j++];

            } else if (right[k].distance > left[j].distance + EPSILON) {
                arr[i++] = right[k++];

            } else {
                if (left[j].index < right[k].index) arr[i++] = left[j++];
                else arr[i++] = right[k++];
            }
        } 
    }

    while (j < leftSize) arr[i++] = left[j++];
    while (k < rightSize) arr[i++] = right[k++];
}

void VectorStore::recursiveMergeSort(Data arr[], int left, int right, bool minimize) const {
    if (left >= right) return;

    int mid = left + (right - left) / 2;
    recursiveMergeSort(arr, left, mid, minimize);
    recursiveMergeSort(arr, mid + 1, right, minimize);
    int leftSize = mid - left + 1;
    int rightSize = right - mid;
    Data* leftArray = new Data[leftSize];
    Data* rightArray = new Data[rightSize];

    for (int i = 0; i < leftSize; ++i) leftArray[i] = arr[left + i];
    for (int j = 0; j < rightSize; ++j) rightArray[j] = arr[mid + 1 + j];
    merge(arr + left, leftArray, leftSize, rightArray, rightSize, minimize);

    delete[] leftArray;
    delete[] rightArray;
}

int* VectorStore::caculatingTopK(const SinglyLinkedList<float> &query, int k, double (VectorStore::*func)(const SinglyLinkedList<float> &, const SinglyLinkedList<float> &) const, bool minimize) const {
    Data* dataMatrix = new Data[count];

    for (int i = 0; i < count; ++i) {
        dataMatrix[i].distance = (this->*func)(query, *(records.get(i)->vector));
        dataMatrix[i].index = i;
    }

    recursiveMergeSort(dataMatrix, 0, count - 1, minimize);

    int* topKIndices = new int[k];
    for (int i = 0; i < k; ++i) {
        topKIndices[i] = dataMatrix[i].index;
    }

    delete[] dataMatrix; 

    return topKIndices;
}

int* VectorStore::topKNearest(const SinglyLinkedList<float> &query, int k, const string &metric) const {
    checkMetric(metric);
    if (k <= 0 || k > count) throw invalid_k_value();

    if (metric == "cosine") return caculatingTopK(query, k, &VectorStore::cosineSimilarity, false);
    else if (metric == "euclidean") return caculatingTopK(query, k, &VectorStore::l2Distance, true);
    else return caculatingTopK(query, k, &VectorStore::l1Distance, true);
}


// Explicit template instantiation for char, string, int, double, float, and Point

template class ArrayList<char>;
template class ArrayList<string>;
template class ArrayList<int>;
template class ArrayList<double>;
template class ArrayList<float>;
template class ArrayList<Point>;
template class ArrayList<VectorStore::VectorRecord*>;

template class SinglyLinkedList<char>;
template class SinglyLinkedList<string>;
template class SinglyLinkedList<int>;
template class SinglyLinkedList<double>;
template class SinglyLinkedList<float>;
template class SinglyLinkedList<Point>;
