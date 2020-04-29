Imports System.Drawing.Imaging

Public Class MDIChlor
    Inherits System.Windows.Forms.Form

#Region " Windows Form Designer generated code "

    Public Sub New()
        MyBase.New()

        'This call is required by the Windows Form Designer.
        InitializeComponent()

        'Add any initialization after the InitializeComponent() call
        frm01 = New frmOption1
    End Sub

    'Form overrides dispose to clean up the component list.
    Protected Overloads Overrides Sub Dispose(ByVal disposing As Boolean)
        If disposing Then
            If Not (components Is Nothing) Then
                components.Dispose()
            End If
        End If
        MyBase.Dispose(disposing)
    End Sub

    'Required by the Windows Form Designer
    Private components As System.ComponentModel.IContainer

    'NOTE: The following procedure is required by the Windows Form Designer
    'It can be modified using the Windows Form Designer.  
    'Do not modify it using the code editor.
    Public WithEvents MainMenu1 As System.Windows.Forms.MainMenu
    Public WithEvents mnuTop As System.Windows.Forms.MenuItem
    Public WithEvents _MnuProject_4 As System.Windows.Forms.MenuItem
    Public WithEvents _MnuProject_5 As System.Windows.Forms.MenuItem
    Public WithEvents ToolTip1 As System.Windows.Forms.ToolTip
    Friend WithEvents MenuItem1 As System.Windows.Forms.MenuItem
    Friend WithEvents MenuItem5 As System.Windows.Forms.MenuItem
    Friend WithEvents MenuItem6 As System.Windows.Forms.MenuItem
    Friend WithEvents MenuItem7 As System.Windows.Forms.MenuItem
    Friend WithEvents MenuItem8 As System.Windows.Forms.MenuItem
    Friend WithEvents MenuItem2 As System.Windows.Forms.MenuItem
    Friend WithEvents MenuItem3 As System.Windows.Forms.MenuItem
    Friend WithEvents MenuItem4 As System.Windows.Forms.MenuItem
    Friend WithEvents MenuItem9 As System.Windows.Forms.MenuItem
    <System.Diagnostics.DebuggerStepThrough()> Private Sub InitializeComponent()
        Me.components = New System.ComponentModel.Container
        Me.MainMenu1 = New System.Windows.Forms.MainMenu
        Me.mnuTop = New System.Windows.Forms.MenuItem
        Me._MnuProject_4 = New System.Windows.Forms.MenuItem
        Me._MnuProject_5 = New System.Windows.Forms.MenuItem
        Me.MenuItem1 = New System.Windows.Forms.MenuItem
        Me.MenuItem2 = New System.Windows.Forms.MenuItem
        Me.MenuItem5 = New System.Windows.Forms.MenuItem
        Me.MenuItem6 = New System.Windows.Forms.MenuItem
        Me.MenuItem7 = New System.Windows.Forms.MenuItem
        Me.MenuItem8 = New System.Windows.Forms.MenuItem
        Me.MenuItem3 = New System.Windows.Forms.MenuItem
        Me.MenuItem4 = New System.Windows.Forms.MenuItem
        Me.ToolTip1 = New System.Windows.Forms.ToolTip(Me.components)
        Me.MenuItem9 = New System.Windows.Forms.MenuItem
        '
        'MainMenu1
        '
        Me.MainMenu1.MenuItems.AddRange(New System.Windows.Forms.MenuItem() {Me.mnuTop, Me.MenuItem1, Me.MenuItem5, Me.MenuItem3})
        '
        'mnuTop
        '
        Me.mnuTop.Index = 0
        Me.mnuTop.MenuItems.AddRange(New System.Windows.Forms.MenuItem() {Me._MnuProject_4, Me._MnuProject_5})
        Me.mnuTop.MergeType = System.Windows.Forms.MenuMerge.Remove
        Me.mnuTop.Text = "Pro&ject"
        '
        '_MnuProject_4
        '
        Me._MnuProject_4.Index = 0
        Me._MnuProject_4.MergeType = System.Windows.Forms.MenuMerge.Remove
        Me._MnuProject_4.Text = "-"
        '
        '_MnuProject_5
        '
        Me._MnuProject_5.Index = 1
        Me._MnuProject_5.MergeType = System.Windows.Forms.MenuMerge.Remove
        Me._MnuProject_5.Text = "E&xit"
        '
        'MenuItem1
        '
        Me.MenuItem1.Index = 1
        Me.MenuItem1.MenuItems.AddRange(New System.Windows.Forms.MenuItem() {Me.MenuItem2})
        Me.MenuItem1.Text = "C&limat"
        '
        'MenuItem2
        '
        Me.MenuItem2.Index = 0
        Me.MenuItem2.Text = "MeteoFile"
        '
        'MenuItem5
        '
        Me.MenuItem5.Index = 2
        Me.MenuItem5.MenuItems.AddRange(New System.Windows.Forms.MenuItem() {Me.MenuItem6, Me.MenuItem7, Me.MenuItem8})
        Me.MenuItem5.Text = "&Transport model"
        '
        'MenuItem6
        '
        Me.MenuItem6.Index = 0
        Me.MenuItem6.Text = "In&put"
        '
        'MenuItem7
        '
        Me.MenuItem7.Index = 1
        Me.MenuItem7.Text = "&Calcul"
        '
        'MenuItem8
        '
        Me.MenuItem8.Index = 2
        Me.MenuItem8.Text = "&Graph"
        '
        'MenuItem3
        '
        Me.MenuItem3.Index = 3
        Me.MenuItem3.MenuItems.AddRange(New System.Windows.Forms.MenuItem() {Me.MenuItem4, Me.MenuItem9})
        Me.MenuItem3.Text = "Probabiliste"
        '
        'MenuItem4
        '
        Me.MenuItem4.Index = 0
        Me.MenuItem4.Text = "Analyse"
        '
        'MenuItem9
        '
        Me.MenuItem9.Index = 1
        Me.MenuItem9.Text = "Graphique"
        '
        'MDIChlor
        '
        Me.AutoScaleBaseSize = New System.Drawing.Size(5, 13)
        Me.BackColor = System.Drawing.SystemColors.AppWorkspace
        Me.ClientSize = New System.Drawing.Size(773, 445)
        Me.FormBorderStyle = System.Windows.Forms.FormBorderStyle.FixedSingle
        Me.IsMdiContainer = True
        Me.Location = New System.Drawing.Point(11, 57)
        Me.Menu = Me.MainMenu1
        Me.Name = "MDIChlor"
        Me.Text = "EPFL  Chloride Ion Penetration Model"
        Me.WindowState = System.Windows.Forms.FormWindowState.Maximized

    End Sub

#End Region

    Public CPTerase As Short = 0
    Dim Para5 As Short
    Dim frm01 As frmOption1
    Dim frm02 As New frmChlor

    'Lorsque la fenêtre est activée
    Private Sub MDIChlor_Load(ByVal eventSender As System.Object, ByVal eventArgs As System.EventArgs) Handles MyBase.Load
        Me.IsMdiContainer = True
        Me.WindowState = FormWindowState.Maximized
        frm02.Left = 0
        frm02.Top = 0
        frm02.Height = (Me.Height)
        frm02.Width = (Me.Width)
        frm02.Hide()
    End Sub

    'quitter le programme
    Private Sub _MnuProject_5_Click(ByVal sender As System.Object, ByVal e As System.EventArgs) Handles _MnuProject_5.Click
        Me.Close()
        End
    End Sub

    'Ouverture de Input dans le menu déroulant
    Private Sub MenuItem6_Click(ByVal sender As System.Object, ByVal e As System.EventArgs) Handles MenuItem6.Click
        frm01.ShowDialog()
        frm01.Hide()
        frm01.Close()
    End Sub

    'Lancement de Open dans le menu déroulant
    Private Sub MenuItem7_Click(ByVal sender As System.Object, ByVal e As System.EventArgs) Handles MenuItem7.Click
        ' Moisture Parameters 
        ' Diffusion Coefficient
        Dim aa As Single ' Function coefficient
        Dim Hc As Single ' Function coefficient
        Dim ab As Single ' Function coefficient
        Dim tc As Single ' Function coefficient
        Dim ImpHydr As Boolean
        ' chloride
        Dim LambdaT(1) As Single
        Dim LambdaH(1) As Single
        Dim aOH As Single
        Dim EbG As Single
        Dim toG As Single
        Dim faG As Single
        Dim hMin As Single
        Dim hEcart As Single
        Dim wMin As Single
        Dim wEcart As Single
        Dim CTmin As Single
        Dim CTecart As Single
        Dim CLmin As Single
        Dim CLecart As Single
        Dim Tecart As Single
        Dim H_snap As Single
        Dim Retard As Single
        ' Structural data
        Dim TimeMax As Single
        Dim Length As Single 'length of the layer [mm]
        Dim Ne As Short 'Number of finite elements
        Dim Le(1) As Decimal 'element length
        Dim PosProf(1) As Decimal
        'Computational data
        Dim Dofs As Short
        'computationals values
        Dim Duration As Single
        Dim DeltaT As Single
        'boundary conditions
        Dim NEXPO As Short
        Dim NQUAL As Short
        Dim taff As Single
        Dim Hsauv As Single
        Dim Wsauv As Single
        Dim CTsauv As Single
        Dim CLsauv As Single
        Dim Tsauv As Single
        Dim Carbsauv As Single
        Dim Filebeton(1) As String
        Dim Fileres(1) As String
        Dim PD(1) As Single
        Dim qGran(1) As Single
        Dim Dcl(1) As Single
        Dim SAT(1) As Single
        Dim ciment(1) As Single
        Dim FileGexpo(1) As String
        Dim FileDexpo(1) As String
        Dim proba(1, 1) As Single
        Dim ChoixRep As Short
        Dim nChmt As Short
        Dim Nbreel(1) As Short
        Dim LenApp(1) As Single
        Dim EC(1) As Single
        Dim tProt(1) As Single
        Dim Vct(1) As Single
        Dim Nct(1) As Single
        Dim capCal As Single
        Dim Hydr(1) As Single
        '---------------------------------------------------------------------------------------
        'Lecture du fichier txt
        Dim Para1 As Single
        Dim Para2 As Single
        Dim Para3 As Single
        Dim Para4 As Single
        Dim test As Single
        Dim i As Short
        Dim j As Short
        Dim nFic As Integer = FreeFile()
        Dim OutFile As String
        Dim Filtre As String
        Dim Index As Short
        Dim Directoire As Boolean
        Dim Titre As String
        Dim Canc As Boolean = False
        Dim ED(1) As Single
        Dim ToHydr(1) As Single
        Dim Ecl(1) As Single
        Dim ToCl(1) As Single
        Dim PostFile As String
        Dim Creadtherm(1, 1) As Double
        Dim Creadhydr(1, 1) As Double
        Dim Creadion(1, 1) As Double
        Dim Ctherm(1) As Double
        Dim Chydr(1) As Double
        Dim Cion(1) As Double
        Dim Nbre1 As Short
        Dim Nbre2 As Short
        Dim Nbre3 As Short
        Dim Var As Double
        Dim Var1 As Double
        Dim GyCO2 As Single
        Dim DyCO2 As Single
        Dim RoC(1) As Single
        Dim RoA(1) As Single

        ''''''''''''''''''''''''''''''''''''''''''''
        Filtre = "Text files (INPUT_*.txt)|INPUT_*.txt"
        Index = 1
        Directoire = True
        Titre = "Sélectionner le fichier d'exposition"
        OpenDialog(OutFile, Canc, Filtre, Index, Directoire, Titre)
        If Canc = True Then GoTo b
        '''''''''''''''''''''''''''''''''''''''''''''

        FileOpen(nFic, OutFile, OpenMode.Input, OpenAccess.Read, OpenShare.Shared)
        FilePost(OutFile, PostFile)

        Input(nFic, Length)
        Input(nFic, Ne)
        Input(nFic, ChoixRep)
        ReDim Le(Ne + CShort(1))
        ReDim PosProf(Ne + CShort(2))

        PosProf(1) = 0
        Select Case ChoixRep
            Case CShort(1)
                For i = CShort(1) To Ne
                    Le(i) = CDec(Length) / CDec(Ne)
                    PosProf(i + CShort(1)) = PosProf(i) + Le(i)
                Next i
            Case CShort(2)
                Input(nFic, Le(1))
                Para1 = CSng(0)
                Para2 = CSng(0)
                For i = CShort(1) To Ne / CShort(2) - CShort(1)
                    Para1 = Para1 + CSng(1)
                    Para2 = Para2 + Para1
                Next i
                If Le(1) > CDec(Length) Then GoTo b
                Para1 = (Length / CSng(2) - CSng(Ne) / CSng(2) * CSng(Le(1))) / Para2
                PosProf(2) = Le(1)
                PosProf(Ne + CShort(1)) = CDec(Length)
                PosProf(Ne) = CDec(Length) - Le(1)
                For i = CShort(2) To Ne / CShort(2)
                    Le(i) = Le(1) + (CDec(i) - CDec(1)) * CDec(Para1)
                    PosProf(i + CShort(1)) = PosProf(i) + Le(i)
                    Le(Ne - i) = Le(i)
                    PosProf(Ne + CShort(1) - i) = PosProf(Ne + CShort(2) - i) - Le(i)
                Next i
            Case 3
                Input(nFic, Le(1))
                Para3 = CSng(Le(1))
                test = CSng(10)
                If Le(1) > CDec(Length) Then GoTo b
                Do While System.Math.Abs(test) > 0.0001
                    Para1 = CSng(1) - Length / (CSng(2) * CSng(Le(1)))
                    Para2 = CSng(1)
                    For i = CShort(1) To Ne / CShort(2) - CShort(1)
                        Para1 = Para1 + Para3 ^ CSng(i)
                        If i < Ne / CShort(2) - CShort(1) Then Para2 = Para2 + (CSng(i) + CSng(1)) * Para3 ^ CSng(i)
                    Next i
                    Para4 = Para3 - Para1 / Para2
                    test = Para4 - Para3
                    Para3 = Para4
                Loop
                PosProf(2) = Le(1)
                PosProf(Ne + CShort(1)) = CDec(Length)
                PosProf(Ne) = CDec(Length) - Le(1)
                For i = CShort(2) To Ne / CShort(2)
                    Le(i) = Le(1) * CDec(Para4) ^ (CDec(i) - CDec(1))
                    PosProf(i + CShort(1)) = PosProf(i) + Le(i)
                    Le(Ne - i) = Le(i)
                    PosProf(Ne + CShort(1) - i) = PosProf(Ne + CShort(2) - i) - Le(i)
                Next i
            Case 4
                Input(nFic, nChmt)
                ReDim Nbreel(nChmt - CShort(1))
                ReDim LenApp(nChmt - CShort(1))
                PosProf(Ne + CShort(1)) = CDec(Length)
                Para1 = CSng(0)
                Para2 = CSng(0)
                For i = CShort(1) To nChmt - CShort(1)
                    Input(nFic, Nbreel(i))
                    Input(nFic, LenApp(i))
                    For j = CShort(1) To Nbreel(i)
                        Le(j + CShort(Para1)) = CDec(LenApp(i)) / CDec(Nbreel(i))
                        PosProf(j + CShort(1) + CShort(Para1)) = PosProf(j + CShort(Para1)) + Le(j + CShort(Para1))
                        Le(Ne + CShort(1) - j - CShort(Para1)) = Le(j + CShort(Para1))
                        PosProf(Ne + CShort(1) - j - CShort(Para1)) = PosProf(Ne + CShort(2) - j - CShort(Para1)) - Le(j + CShort(Para1))
                    Next j
                    Para1 = Para1 + CSng(Nbreel(i))
                    Para2 = Para2 + LenApp(i)
                Next i
                Para3 = CSng(Ne) - CSng(2) * Para1
                Para2 = Length - CSng(2) * Para2
                For i = CShort(Para1) + CShort(1) To CShort(Para1) + CShort(Para3)
                    Le(i) = CDec(Para2) / CDec(Para3)
                    PosProf(i + CShort(1)) = PosProf(i) + Le(i)
                Next i
            Case 5
                Input(nFic, nChmt)
                ReDim Nbreel(nChmt - CShort(1))
                ReDim LenApp(nChmt - CShort(1))
                PosProf(Ne + CShort(1)) = CDec(Length)
                Para1 = CSng(0)
                Para2 = CSng(0)
                For i = CShort(1) To nChmt - CShort(1)
                    Input(nFic, Nbreel(i))
                    Input(nFic, LenApp(i))
                    For j = CShort(1) To Nbreel(i)
                        Le(j + CShort(Para1)) = CDec(LenApp(i)) / CDec(Nbreel(i))
                        PosProf(j + CShort(1) + CShort(Para1)) = PosProf(j + CShort(Para1)) + Le(j + CShort(Para1))
                    Next j
                    Para1 = Para1 + CSng(Nbreel(i))
                    Para2 = Para2 + LenApp(i)
                Next i
                Para3 = CSng(Ne) - Para1
                Para2 = Length - Para2
                For i = CShort(Para1) + CShort(1) To CShort(Para1) + CShort(Para3)
                    Le(i) = CDec(Para2) / CDec(Para3)
                    PosProf(i + CShort(1)) = PosProf(i) + Le(i)
                Next i
        End Select

        Input(nFic, TimeMax)
        Input(nFic, DeltaT)
        Input(nFic, taff)
        Input(nFic, Hsauv)
        Input(nFic, Wsauv)
        Input(nFic, CTsauv)
        Input(nFic, CLsauv)
        Input(nFic, Tsauv)
        Input(nFic, Carbsauv)
        Input(nFic, hMin)
        Input(nFic, hEcart)
        Input(nFic, wMin)
        Input(nFic, wEcart)
        Input(nFic, CLmin)
        Input(nFic, CLecart)
        Input(nFic, CTmin)
        Input(nFic, CTecart)
        Input(nFic, Tecart)
        Input(nFic, aa)
        Input(nFic, Hc)
        Input(nFic, ab)
        Input(nFic, tc)
        Input(nFic, ImpHydr)
        Input(nFic, H_snap)
        Input(nFic, Retard)
        Input(nFic, aOH)
        Input(nFic, EbG)
        Input(nFic, toG)
        Input(nFic, faG)
        Input(nFic, capCal)
        Input(nFic, GyCO2)
        Input(nFic, DyCO2)

        Input(nFic, NEXPO)
        ReDim FileGexpo(NEXPO)
        ReDim FileDexpo(NEXPO)
        For i = CShort(1) To NEXPO
            Input(nFic, FileGexpo(i))
            Input(nFic, FileDexpo(i))
        Next i

        Input(nFic, NQUAL)
        ReDim Filebeton(NQUAL)
        ReDim Fileres(NQUAL)
        ReDim PD(NQUAL)
        ReDim Dcl(NQUAL)
        ReDim qGran(NQUAL)
        ReDim LambdaH(NQUAL)
        ReDim LambdaT(NQUAL)
        ReDim SAT(NQUAL)
        ReDim ciment(NQUAL)
        ReDim EC(NQUAL)
        ReDim tProt(NQUAL)
        ReDim Vct(NQUAL)
        ReDim Nct(NQUAL)
        ReDim Hydr(NQUAL)
        ReDim ED(NQUAL)
        ReDim ToHydr(NQUAL)
        ReDim Ecl(NQUAL)
        ReDim ToCl(NQUAL)
        ReDim RoA(NQUAL)
        ReDim RoC(NQUAL)
        ReDim proba(NQUAL, 19)
        For i = CShort(1) To NQUAL
            Input(nFic, Filebeton(i))
            Input(nFic, Fileres(i))
            Input(nFic, PD(i))
            Input(nFic, Dcl(i))
            Input(nFic, qGran(i))
            Input(nFic, LambdaH(i))
            Input(nFic, LambdaT(i))
            Input(nFic, SAT(i))
            Input(nFic, ciment(i))
            Input(nFic, EC(i))
            Input(nFic, tProt(i))
            Input(nFic, Hydr(i))
            Input(nFic, Vct(i))
            Input(nFic, Nct(i))
            Input(nFic, ED(i))
            Input(nFic, ToHydr(i))
            Input(nFic, Ecl(i))
            Input(nFic, ToCl(i))
            Input(nFic, RoA(i))
            Input(nFic, RoC(i))
            For j = 0 To 19
                Input(nFic, proba(i, j))
            Next
        Next i

        Input(nFic, Nbre1)
        ReDim Creadtherm(1, Nbre1)
        For i = 1 To Nbre1
            Input(nFic, Creadtherm(0, i))
            Input(nFic, Creadtherm(1, i))
        Next

        Input(nFic, Nbre2)
        ReDim Creadhydr(1, Nbre2)
        For i = 1 To Nbre2
            Input(nFic, Creadhydr(0, i))
            Input(nFic, Creadhydr(1, i))
            Creadhydr(1, i) = Creadhydr(1, i) / 100
        Next

        Input(nFic, Nbre3)
        ReDim Creadion(1, Nbre3)
        For i = 1 To Nbre3
            Input(nFic, Creadion(0, i))
            Input(nFic, Creadion(1, i))
        Next

        Dofs = Ne + CShort(1)
        'calcul des conditions initiales
        j = 0
        ReDim Ctherm(Dofs + 1)
        For i = 1 To Nbre1 - 1
            Var = (Creadtherm(1, i + 1) - Creadtherm(1, i)) / (Creadtherm(0, i + 1) - Creadtherm(0, i))
            Var1 = Creadtherm(1, i) - Var * Creadtherm(0, i)
            Do While PosProf(j) <= Creadtherm(0, i + 1)
                Ctherm(j) = Var * PosProf(j) + Var1
                j = j + 1
                If j > Dofs Then Exit Do
            Loop
            Ctherm(Dofs + 1) = Ctherm(Dofs)
        Next
        j = 0
        ReDim Chydr(Dofs + 1)
        For i = 1 To Nbre2 - 1
            Var = (Creadhydr(1, i + 1) - Creadhydr(1, i)) / (Creadhydr(0, i + 1) - Creadhydr(0, i))
            Var1 = Creadhydr(1, i) - Var * Creadhydr(0, i)
            Do While PosProf(j) <= Creadhydr(0, i + 1)
                Chydr(j) = Var * PosProf(j) + Var1
                j = j + 1
                If j > Dofs Then Exit Do
            Loop
            Chydr(Dofs + 1) = Chydr(Dofs)
            If j > Dofs + 1 Then Exit For
        Next
        j = 0
        ReDim Cion(Dofs + 1)
        For i = 1 To Nbre3 - 1
            Var = (Creadion(1, i + 1) - Creadion(1, i)) / (Creadion(0, i + 1) - Creadion(0, i))
            Var1 = Creadion(1, i) - Var * Creadion(0, i)
            Do While PosProf(j) <= Creadion(0, i + 1)
                Cion(j) = Var * PosProf(j) + Var1
                j = j + 1
                If j > Dofs Then Exit Do
            Loop
            Cion(Dofs + 1) = Cion(Dofs)
        Next

        Le(0) = CDec(1)       'couche limite
        Le(Dofs) = CDec(1)    'couche limite

        FileClose(nFic)
        '---------------------------------------------------------------------------------------
        Dim myThread As System.Threading.Thread

        SetParameters(frm02, Length, Ne, ChoixRep, Le, PosProf, nChmt, Nbreel, LenApp, TimeMax, DeltaT, taff, Hsauv, Wsauv, CTsauv, CLsauv, Tsauv, Carbsauv, hMin, hEcart, wMin, wEcart, CLmin, CLecart, CTmin, CTecart, Tecart, aa, Hc, ab, tc, ImpHydr, H_snap, Retard, aOH, EbG, toG, faG, NEXPO, FileGexpo, FileDexpo, NQUAL, Filebeton, Fileres, PD, Dcl, capCal, LambdaH, LambdaT, SAT, ciment, EC, tProt, Vct, Nct, proba, Dofs, qGran, Hydr, ED, ToHydr, Ecl, ToCl, PostFile, Ctherm, Chydr, Cion, GyCO2, DyCO2, RoA, RoC)

        myThread = New System.Threading.Thread(AddressOf Compute_All)

        frm02.MdiParent = Me
        frm02.Show()
        If Para5 <> CShort(1) Then
            myThread.Start()
        End If

b:      'user pressed cancel error

    End Sub

    'construction de graphique
    Private Sub MenuItem8_Click(ByVal sender As System.Object, ByVal e As System.EventArgs) Handles MenuItem8.Click
        Dim frm03 As frmScale1

        frm03 = New frmScale1
        frm02.MdiParent = Me
        frm02.Show()
        manage_graph(frm02, frm03)
    End Sub

    'traitement des données météorologiques
    Private Sub MenuItem2_Click(ByVal sender As System.Object, ByVal e As System.EventArgs) Handles MenuItem2.Click
        MeteoTreatment()
    End Sub

    'lancement de l'approche probabiliste
    Private Sub MenuItem4_Click(ByVal sender As System.Object, ByVal e As System.EventArgs) Handles MenuItem4.Click
        Dim frm03 As frmProb

        frm03 = New frmProb
        frm03.MdiParent = Me
        frm03.Show()

    End Sub

    'lancement du traitement graphique des probabilités
    Private Sub MenuItem9_Click(ByVal sender As System.Object, ByVal e As System.EventArgs) Handles MenuItem9.Click
        Dim frm03 As frmScale1

        frm03 = New frmScale1
        frm02.MdiParent = Me
        frm02.Show()
        ProbGraph(frm02, frm03)

    End Sub

End Class
