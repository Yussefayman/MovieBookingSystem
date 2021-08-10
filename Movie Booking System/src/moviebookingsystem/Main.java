package moviebookingsystem;
	
import javafx.application.Application;
import javafx.fxml.FXMLLoader;
import javafx.stage.Stage;
import javafx.scene.Parent;
import javafx.scene.Scene;
import javafx.scene.layout.BorderPane;


public class Main extends Application {
  static Stage window ;
	@Override
	 public void start(Stage stage) throws Exception {
             window=stage;
        Parent root = FXMLLoader.load(getClass().getResource("LoginGui.fxml"));

        Scene scene = new Scene(root);

        stage.setScene(scene);
        stage.setTitle("Registration Form");
        stage.show();
    }

	public static void main(String[] args) {
		launch(args);
	}
}
