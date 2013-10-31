#ifndef MODIFIER_H
#define MODIFIER_H

class MoleculeSystem;

class Modifier
{
public:
    Modifier(MoleculeSystem* moleculeSystem);
    virtual ~Modifier() {}

    virtual void apply() = 0;

protected:
    MoleculeSystem* m_moleculeSystem;
};

#endif // MODIFIER_H
